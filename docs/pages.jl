# This file assumes `dir` is the directory for the package! dir = @__DIR__() * "/.."

dir = @__DIR__() * "/.."

cp(joinpath(dir, "markdown"), joinpath(dir, "docs", "src"), force=true)
cp(joinpath(dir, "docs", "extrasrc", "assets"), joinpath(dir, "docs", "src", "assets"), force=true)
cp(joinpath(dir, "README.md"), joinpath(dir, "docs", "src", "index.md"), force=true)
tutorialsdir = joinpath(dir, "docs", "src")

pages = Any["SciMLTutorials.jl: Tutorials for Scientific Machine Learning (SciML), Equation Solvers, and AI for Science"=>"index.md"]

for folder in readdir(tutorialsdir)
    newpages = Any[]
    if folder[end-2:end] != ".md" && folder != "Testing" && folder != "figures" && folder != "assets"
        for file in filter(x -> x[end-2:end] == ".md", readdir(
            joinpath(tutorialsdir, folder)))
            try
                filecontents = readlines(joinpath(tutorialsdir, folder, file))
                title = filecontents[3][9:end-1]

                # Cut out the first 5 lines from the file to remove the Weave header stuff
                open(joinpath(tutorialsdir, folder, file), "w") do output
                    println(output, "# $title")
                    for line in Iterators.drop(filecontents, 4)
                        println(output, line)
                    end
                end
                push!(newpages, title => joinpath(folder, file))
            catch e
                @show folder, file, e
            end
        end
        push!(pages, folder => newpages)
    end
end

# The result is in alphabetical order, change to the wanted order

permute!(pages,
    [1, 4, 6, 2, 3, 5]
)

names = [
    "SciMLTutorials.jl: Tutorials for Scientific Machine Learning (SciML) and Equation Solvers",
    "Ordinary Differential Equation (ODE) Examples",
    "Special Analyses of ODEs",
    "Advanced ODE Examples",
    "Symbolic-Numeric Approaches",
    "Workshop Exercises"]

for i in 1:length(pages)
    pages[i] = names[i] => pages[i][2]
end
