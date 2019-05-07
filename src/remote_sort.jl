"""
Transfers all data of type `T` under the current directory to `remote`
"""
function transfer_data(::Type{T};remote="nus-hpc:/hpctmp/lsihr/", dryrun=false, start_sort=true) where T <: DPHT.DPHData
    fname = DPHT.filename(T)
    ll = DPHT.level(T)
    lt = DPHT.level()
    if ll == lt
        dirs = ["."]
    else
        dirs = DPHT.get_level_dirs(DPHT.level(T))
    end
    rpath = DPHT.get_relative_path("subjects")
    files = [joinpath(".", d,fname) for d in dirs]
    if dryrun
        ss = ["--dry-run"]
    else
        ss = []
    end
    hostname,pth = split(remote, ":")
    pth = joinpath(pth, rpath)
    run(`ssh $hostname mkdir -p $pth`)
    for f in files
        cmd =`rsync -ruvLR --partial --progress $ss $f $(joinpath(remote, rpath,""))`
        if dryrun
            @show files
            @show cmd
            @show pth
        end
        run(cmd)
        _pth,_fname = splitdir(f)
        _pth = joinpath(pth, _pth)
        if start_sort
            run(`ssh $hostname cd $_pth; ~/programming/hmmsort/hmmsort_pbs.py ~/programming/hmmsort`)
        end
    end
end
