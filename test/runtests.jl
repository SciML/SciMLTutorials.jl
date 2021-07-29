using SciMLTutorials
tutorials_dir = joinpath(dirname(@__DIR__), "tutorials")
SciMLTutorials.weave_file(joinpath(tutorials_dir, "Testing"), "test.jmd")
