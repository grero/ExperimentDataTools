installed_packages = keys(Pkg.installed())
isinstalled(pkg) = pkg in installed_packages

if !isinstalled("SpikeExtraction")
    Pkg.clone("https://github.com/grero/SpikeExtraction.jl.git", "SpikeExtraction")
end
if !isinstalled("HMMSpikeSorter")
    Pkg.clone("https://github.com/grero/HMMSpikeSorter.jl.git", "HMMSpikeSorter")
end
if !isinstalled("Eyelink")
    Pkg.clone("https://github.com/grero/Eyelink.jl.git", "Eyelink")
end
if !isinstalled("Stimulus")
    Pkg.clone("https://bitbucket.org/rherikstad/stmulus.jl.git", "Stimulus")
end
if !isinstalled("Spiketrains")
    Pkg.clone("https://bitbucket.org/rherikstad/spiketrains.jl.git", "Spiketrains")
end
if !isinstalled("LFPTools")
    Pkg.clone("https://github.com/grero/LFPTools.jl.git", "LFPTools")
end
if !isinstalled("RippleTools")
    Pkg.clone("http://gitlab.thornbush/grero/RippleTools.jl.git", "RippleTools")
end
if !isinstalled("LegacyFileReaders")
    Pkg.clone("https://github.com/grero/LegacyFileReaders.jl.git", "LegacyFileReaders")
end
if !isinstalled("DataProcessingHierarchyTools")
    Pkg.clone("https://github.com/grero/DataProcessingHierarchyTools.jl.git", "DataProcessingHierarchyTools")
end
Pkg.update()
if isinstalled("ExperimentDataTools")
    Pg.rm("ExperimentDataTools")
end
Pkg.clone(pwd())
Pkg.build("ExperimentDataTools")
