# Introduction
Various tools to process rawdata into lowpass,highpass and multiunit streams.

## Creating sessions

To copy data for the various sessions into session sub-directories

```julia
;ls
w5_16.nev            w5_16_2.edf          w5_16_3_settings.txt
w5_16.ns5            w5_16_2_results.txt  w5_16_4.edf
w5_16_1.edf          w5_16_2_settings.txt w5_16_4_results.txt
w5_16_1_results.txt  w5_16_3.edf          w5_16_4_settings.txt
w5_16_1_settings.txt w5_16_3_results.txt

ExperimentDataTools.process()
;ls
session01 session02 session03 session04 w5_16.nev w5_16.ns5
```

## Extract trial markers

To extract trial markers from the event data and those as event_markers.csv

```julia
trials = ExperimentDataTools.Trials()
```

## Processing raw data

To extract lowpass and highpass signals from raw data
```julia
ExperimentDataTools.process_rawdata(File(format"NSHR", "w5_16.ns5"))
```

## Extract multi-unit activity
