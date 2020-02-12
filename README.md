# Chromatin SE analysis pipeline.
- Currently designed as ChipSeq



## To Do before  deployment:
- Toil is unable to locate user-installed programs/ packages such as
  - **R package spp** [for PEAKSQC step (NSR, RSC, Phantom Quality calculation)]. __https://github.com/crazyhottommy/phantompeakqualtools.git__
  - **sicer** [SICERv1.1 with hardcoded path is still utilized, but will prefer SICER2 for performance and support]. __https://github.com/zanglab/SICER2.git__
