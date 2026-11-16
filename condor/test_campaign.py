#!/usr/bin/env python3
"""Standalone tests for condor/campaign.py (no pytest, no Condor and no EOS needed)

Run from the repository root:

    python condor/test_campaign.py

Covers loading and validating the campaign configs, the condor_submit,
make_job_args.py and merge command lines built from them, the files-per-job cap,
the refusal of collision data unless a config allows it, which samples a full or a
partial merge takes, and the bookkeeping of retry rounds and of input files that
keep failing, with synthetic fixtures.
"""

import copy
import os
import re
import sys
import tempfile

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import campaign  # noqa: E402  pylint: disable=wrong-import-position

_PASS = 0
_FAIL = 0


def check(name, cond):
    """Count a passed or failed check, printing the failures"""
    global _PASS, _FAIL
    if cond:
        _PASS += 1
    else:
        _FAIL += 1
        print("FAIL:", name)


def raises(cfg):
    """True if building a Campaign from cfg raises ValueError"""
    try:
        campaign.Campaign(cfg)
    except ValueError:
        return True
    return False


def corrections_config():
    """The corrections_v1 config as a dict"""
    with open(os.path.join(campaign.CONFIG_DIR, "corrections_v1.yaml"), encoding="utf8") as f:
        return yaml.safe_load(f)


def test_configs_load():
    """Both committed configs load; the corrections one has the values that ran"""
    os.environ.setdefault("USER", "someone")
    camp = campaign.load_campaign("corrections_v1")
    user = os.environ["USER"]
    check("name", camp.name == "corrections_v1")
    check("$USER expanded", camp.chunk_dir == f"/store/user/{user}/sidm_condor/CorrectionsChunks_v1")
    check("chunk URL", camp.chunk_url == f"root://cmseos.fnal.gov//store/user/{user}/sidm_condor/CorrectionsChunks_v1")
    check("corrections list", camp.corrections == "pileup,l1_prefiring,muon_id,electron_reco")
    check("monitored", camp.monitor == ["photon_pixel_veto"])
    check("no cap for the record", camp.max_files_per_job is None)
    check("three job lists", len(camp.job_list_specs) == 3)
    check("job lists exist", all(os.path.exists(p) for p in camp.job_lists()[:3]))
    check("logs dir", camp.logs_dir == os.path.join(campaign.CONDOR, "logs_corrections"))
    check("by path too", campaign.load_campaign(os.path.join(campaign.CONFIG_DIR, "corrections_v1.yaml")).name
          == "corrections_v1")
    template = campaign.load_campaign("template")
    check("template is valid", template.max_files_per_job == 100 and template.corrections == "default")


def test_validation():
    """Typos and inconsistent options are refused"""
    good = corrections_config()
    check("good config builds", not raises(good))
    for key in ("name", "chunk_dir", "channels", "job_lists"):
        cfg = copy.deepcopy(good)
        del cfg[key]
        check(f"missing {key}", raises(cfg))
    cfg = copy.deepcopy(good)
    cfg["chanels"] = cfg.pop("channels")
    check("misspelt key", raises(cfg))
    cfg = copy.deepcopy(good)
    cfg["hist_mode"] = "weighed"
    check("unknown hist_mode", raises(cfg))
    cfg = copy.deepcopy(good)
    cfg["corrections"] = "none"
    check("syst without corrections", raises(cfg))
    cfg = copy.deepcopy(good)
    cfg["chunk_dir"] = "/eos/uscms/store/user/x"
    check("chunk_dir outside /store", raises(cfg))
    cfg = copy.deepcopy(good)
    cfg["merged_dir"] = "/store/user/$SIDM_UNSET_VARIABLE_FOR_TEST/x"
    check("unexpanded variable", raises(cfg))
    cfg = copy.deepcopy(good)
    cfg["max_files_per_job"] = 0
    check("zero cap", raises(cfg))
    cfg = copy.deepcopy(good)
    del cfg["job_lists"][0]["dataset"]
    check("job list without dataset", raises(cfg))
    cfg = copy.deepcopy(good)
    cfg["job_lists"][1]["job_args"] = cfg["job_lists"][0]["job_args"]
    check("two lists in one file", raises(cfg))
    cfg = copy.deepcopy(good)
    cfg["chunk_dir"] = "root://cmseos.fnal.gov//store/user/x/chunks"
    check("redirector stripped", campaign.Campaign(cfg).chunk_dir == "/store/user/x/chunks")
    cfg = copy.deepcopy(good)
    cfg["channels"] = "base, 2mu2e"
    check("comma-separated names", campaign.Campaign(cfg).channels == ["base", "2mu2e"])


def test_submit():
    """The submit macros reproduce the arguments the corrections campaign ran with"""
    camp = campaign.load_campaign("corrections_v1")
    job_list = os.path.join(campaign.CONDOR, "job_args_corrections_2mu2e.txt")
    macros = camp.submit_macros(job_list)
    with open(os.path.join(campaign.CONDOR, campaign.SUBMIT_FILE), encoding="utf8") as f:
        sub = "\n".join(line for line in f if not line.lstrip().startswith("#"))
    used = set(re.findall(r"\$\((\w+)\)", sub)) - {"sample", "chunk", "filelist", "Cluster", "Process"}
    check("every macro of the submit file is given", used == set(macros))
    arguments = re.search(r"^arguments = (.*)$", sub, re.M).group(1)
    expanded = re.sub(r"\$\((\w+)\)", lambda m: macros.get(m.group(1), m.group(0)), arguments)
    # the nine arguments of the old submit_corrections.sub, then the selections file and the data guard
    old = ("$(sample) $(chunk) $(filelist) root://cmseos.fnal.gov//store/user/$ENV(USER)/sidm_condor/"
           "CorrectionsChunks_v1 base,2mu2e,4mu,2mu2e_sr,4mu_sr,Z_to_mumu_control_region "
           "corrections_study pileup,l1_prefiring,muon_id,electron_reco syst photon_pixel_veto")
    check("arguments as in submit_corrections.sub", expanded == old.replace("$ENV(USER)", os.environ["USER"])
          + " configs/selections.yaml no")
    template = campaign.load_campaign("template").submit_macros(job_list)
    check("argument 9 is none, not empty", template["MONITOR"] == "none")
    check("data guard kept", template["ALLOW_DATA"] == "no")
    with open(os.path.join(campaign.CONDOR, "run_job.sh"), encoding="utf8") as f:
        run_job = f.read()
    check("run_job.sh takes the 11th argument", "${11:-no}" in run_job and "--allow-data" in run_job)
    check("job list relative to condor/", macros["JOB_ARGS"] == "job_args_corrections_2mu2e.txt")
    check("logs", macros["LOGS"] == "logs_corrections")
    check("memory", macros["MEMORY"] == "2GB" and camp.submit_macros(job_list, "8GB")["MEMORY"] == "8GB")
    check("stall guard", macros["STALL_HOURS"] == "12" and macros["STALL_RESTARTS"] == "1")
    cmd = camp.submit_command(job_list)
    check("command ends with the submit file", cmd[-1] == campaign.SUBMIT_FILE and "CORRECTIONS=" + camp.corrections in cmd)


def test_values():
    """The selections file, its channels, and values that would break the arguments line"""
    good = corrections_config()
    for key, value in (("selections_cfg", "configs/no_such_file.yaml"),
                       ("selections_cfg", os.path.join(campaign.REPO, "sidm", "configs", "selections.yaml")),
                       ("channels", ["base", "no_such_channel"]),
                       ("channels", ["base", "2mu2e sr"]),
                       ("hist_collections", ["corrections_study'"]),
                       ("chunk_dir", "/store/user/x/my chunks"),
                       ("request_memory", "2 GB"),
                       ("allow_data", "yes")):
        cfg = copy.deepcopy(good)
        cfg[key] = value
        check(f"{key} = {value!r} refused", raises(cfg))
    cfg = copy.deepcopy(good)
    cfg["selections_cfg"] = "configs/selections.yaml"
    check("selections file given explicitly", campaign.Campaign(cfg).selections_cfg == "configs/selections.yaml")


def test_data_campaigns():
    """A campaign over collision data needs allow_data: true"""
    saved = campaign.NTUPLE_DIR
    with tempfile.TemporaryDirectory() as tmp:
        try:
            campaign.NTUPLE_DIR = tmp
            menus = {"data_fake.yaml": {"v1": {"samples": {"DoubleMuon_A": {}}}},
                     "mc_fake.yaml": {"v1": {"samples": {"S": {"is_data": False}}}},
                     "flagged_fake.yaml": {"v1": {"samples": {"D": {"is_data": True}}}}}
            for name, menu in menus.items():
                with open(os.path.join(tmp, name), "w", encoding="utf8") as f:
                    yaml.safe_dump(menu, f)
            cfg = corrections_config()
            cfg["job_lists"] = [{"job_args": f"job_args_{name}.txt", "samples_file": "x.txt", "dataset": "v1",
                                 "location_cfg": name} for name in menus]
            camp = campaign.Campaign(cfg)
            check("data found by name and by is_data",
                  camp.data_job_lists() == ["job_args_data_fake.yaml.txt", "job_args_flagged_fake.yaml.txt"])
            try:
                camp.check_data()
                check("data refused without allow_data", False)
            except SystemExit:
                check("data refused without allow_data", True)
            cfg["allow_data"] = True
            allowed = campaign.Campaign(cfg)
            allowed.check_data()
            check("allow_data passes --allow-data", allowed.submit_macros(os.path.join(tmp, "j.txt"))["ALLOW_DATA"]
                  == "allow-data")
            cfg["allow_data"] = False
            cfg["job_lists"] = cfg["job_lists"][1:2]
            mc = campaign.Campaign(cfg)
            mc.check_data()
            check("simulation needs nothing", not mc.data_job_lists())
        finally:
            campaign.NTUPLE_DIR = saved


def test_merge():
    """The merge command, and which samples a full or a partial merge takes"""
    camp = campaign.load_campaign("corrections_v1")
    cmd = camp.merge_command("samples.txt", "filelists", "/store/user/x/merged")
    check("selections file to the sidecar", cmd[cmd.index("--selections-cfg") + 1] == "configs/selections.yaml")
    check("weighted campaign", "--unweighted-hist" not in cmd)
    template = campaign.load_campaign("template")
    check("unweighted campaign", "--unweighted-hist" in template.merge_command("s", "f", "/store/user/x/m"))

    def result(sample, chunk, state):
        return campaign.cl.ChunkResult(sample=sample, chunk=chunk, key=f"{sample}_{chunk}", state=state)
    cl = campaign.cl
    results = [result("A", 0, cl.DONE), result("A", 1, cl.DONE),
               result("B", 0, cl.DONE), result("B", 1, cl.IN_PROGRESS),
               result("C", 0, cl.ZERO_BYTE_OUTPUT), result("C", 1, cl.DONE),
               result("D", 0, cl.IN_PROGRESS), result("D", 1, cl.OOM)]
    chosen, done, left_out = campaign.merge_plan(results)
    check("full merge: complete samples", chosen == ["A"] and done == {"A_0", "A_1"} and not left_out)
    check("full merge: not merged twice", campaign.merge_plan(results, merged={"A"})[0] == [])
    for wanted in (["B"], ["X"]):
        try:
            campaign.merge_plan(results, wanted)
            check(f"full merge refuses {wanted}", False)
        except ValueError:
            check(f"full merge refuses {wanted}", True)
    chosen, done, left_out = campaign.merge_plan(results, partial=True)
    check("partial merge: every sample with a done chunk", chosen == ["A", "B"])
    check("partial merge: only the done chunks", done == {"A_0", "A_1", "B_0"})
    check("partial merge: zero-byte and empty samples left out", sorted(left_out) == ["C", "D"])
    check("partial merge of one sample", campaign.merge_plan(results, ["B"], partial=True)[:2] == (["B"], {"B_0"}))


def test_make_job_args():
    """make_job_args.py command lines and the files-per-job cap"""
    camp = campaign.load_campaign("corrections_v1")
    bkg = camp.make_job_args_command(camp.job_list_specs[2])
    check("map passed", bkg[bkg.index("--files-per-job-map") + 1] == "condor/files_per_job_corrections.yaml")
    check("census skip passed", bkg[bkg.index("--census-skip") + 1] == "backgrounds_skimmed.skip.json")
    check("xcache replaced", "--replace-xcache" in bkg)
    check("no cap for the record", "--max-files-per-job" not in bkg)
    check("filelists dir", bkg[bkg.index("--outdir") + 1] == "condor/filelists_corrections")
    cfg = corrections_config()
    del cfg["max_files_per_job"]
    capped = campaign.Campaign(cfg).make_job_args_command(camp.job_list_specs[2])
    check("default cap of 100", capped[capped.index("--max-files-per-job") + 1] == "100")
    import make_job_args  # pylint: disable=import-outside-toplevel
    fpj = {"tiny": 300, "big": 6}
    check("map value", make_job_args.files_per_job_for("big", fpj, 25) == 6)
    check("default", make_job_args.files_per_job_for("signal", fpj, 25) == 25)
    check("capped", make_job_args.files_per_job_for("tiny", fpj, 25, cap=100) == 100)
    check("below the cap", make_job_args.files_per_job_for("big", fpj, 25, cap=100) == 6)


def test_bookkeeping():
    """Retry rounds, files that failed every attempt, and the bad-file lists"""
    camp = campaign.load_campaign("corrections_v1")
    saved = campaign.CONDOR
    with tempfile.TemporaryDirectory() as tmp:
        try:
            campaign.CONDOR = tmp
            check("first round", campaign.next_round(camp) == 1)
            for name in ("job_args_corrections_retry1.txt", "job_args_corrections_retry2_oom.txt",
                         "job_args_corrections_extra_round3.txt", "job_args_other_retry9.txt"):
                with open(os.path.join(tmp, name), "w", encoding="utf8"):
                    pass
            check("next round after retries and extra rounds", campaign.next_round(camp) == 4)
        finally:
            campaign.CONDOR = saved

        bad_file = "root://cmseos.fnal.gov//store/group/x/bad.root"
        err = os.path.join(tmp, "job.err")
        with open(err, "w", encoding="utf8") as f:
            f.write(f"Failed processing file: WorkItem(dataset='S', filename='{bad_file}', treename='Events')\n" * 2)

        def result(n_failures):
            log = os.path.join(tmp, f"job{n_failures}.log")
            with open(log, "w", encoding="utf8") as f:
                f.write("005 (1.0.0) Job terminated.\n\t(1) Normal termination (return value 1)\n" * n_failures)
            return campaign.cl.ChunkResult(sample="S", chunk=0, key="S_0", state=campaign.cl.UNKNOWN,
                                           log_paths={"log": log, "err": err})
        check("failed every attempt: dropped", campaign.failing_file(result(2)) == bad_file)
        check("failed once: kept", campaign.failing_file(result(1)) is None)
        with open(err, "a", encoding="utf8") as f:
            f.write("Failed processing file: WorkItem(dataset='S', filename='root://h//store/other.root')\n")
        check("two different files: kept", campaign.failing_file(result(2)) is None)

        filelist = os.path.join(tmp, "S_0.txt")
        with open(filelist, "w", encoding="utf8") as f:
            f.write(f"{bad_file}\nroot://cmseos.fnal.gov//store/group/x/good.root\n")
        bad_list = os.path.join(tmp, "bad_files.txt")
        with open(bad_list, "w", encoding="utf8") as f:
            f.write(f"{bad_file}  # S_0: failed every attempt\n")
        bad = campaign.read_file_list(bad_list)
        check("bad-file list read as /store paths", bad == {"/store/group/x/bad.root"})
        check("filelist keeps the good file", campaign.drop_files(filelist, bad)
              and campaign.cl.read_filelist(filelist) == ["root://cmseos.fnal.gov//store/group/x/good.root"])
        check("nothing left", not campaign.drop_files(filelist, {"/store/group/x/good.root"}))
        check("no list, no files", campaign.read_file_list(os.path.join(tmp, "absent.txt")) == set())


if __name__ == "__main__":
    test_configs_load()
    test_validation()
    test_values()
    test_submit()
    test_data_campaigns()
    test_merge()
    test_make_job_args()
    test_bookkeeping()
    print(f"{_PASS} passed, {_FAIL} failed")
    sys.exit(1 if _FAIL else 0)
