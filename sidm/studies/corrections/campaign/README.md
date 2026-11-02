# corrections_v1 campaign records

Bookkeeping of the Condor campaign described in `../README.md` ("The Condor campaign").

* `bad_files.txt`: the 13 input files dropped because they are corrupt (they failed in
  every attempt, and one to four branches each do not decompress), with the chunk each was in.
* `empty_files.txt`: input files with no events.
* `local_chunks.txt`: the four chunks that were processed on an interactive node, with the
  code of the job tarball and the same options, because their Condor jobs kept landing
  on slow worker nodes.

The job lists are `condor/job_args_corrections_*.txt` and the submit file
`condor/submit_corrections.sub`; `sidm/scripts/corrections_campaign.py` keeps its working
state (logs, reports, retry lists) under `condor/campaigns/corrections_v1/`, which is not
committed.
