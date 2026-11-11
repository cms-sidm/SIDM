# corrections_v1 campaign records

Bookkeeping of the Condor campaign described in `../README.md` ("The Condor campaign").

* `bad_files.txt`: the 13 input files dropped because they are corrupt (they failed in
  every attempt, and one to four branches each do not decompress), with the chunk each was in.
* `empty_files.txt`: input files with no events.
* `local_chunks.txt`: the four chunks that were processed on an interactive node, with the
  code of the job tarball and the same options, because their Condor jobs kept landing
  on slow worker nodes.

The campaign is `condor/campaign_configs/corrections_v1.yaml`, with the job lists
`condor/job_args_corrections_*.txt`. `condor/campaign.py` keeps its working state
(bad-file lists, the files skipped per chunk, merge scratch space) under
`condor/campaigns/corrections_v1/`, which is not committed.
