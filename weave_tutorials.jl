using SciMLTutorials
target = ARGS[1]
if isdir(target)
    if !isfile(joinpath(target, "Project.toml"))
        error("Cannot weave folder $(target) without Project.toml!")
    end
    println("Weaving the $(target) folder")
    SciMLTutorials.weave_folder(target)
elseif isfile(target)
    folder = dirname(target)
    file = basename(target)
    println("Weaving $(folder)/$(file)")
    SciMLTutorials.weave_file(folder, file)
else
    error("Unable to find weaving target $(target)!")
end
