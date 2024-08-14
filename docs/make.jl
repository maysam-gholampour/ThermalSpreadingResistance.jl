using ThermalSpreadingResistance
using Documenter

DocMeta.setdocmeta!(ThermalSpreadingResistance, :DocTestSetup, :(using ThermalSpreadingResistance); recursive=true)

makedocs(;
    modules=[ThermalSpreadingResistance],
    authors="maysam-gholampour <meysam.gholampoor@gmail.com> and contributors",
    sitename="ThermalSpreadingResistance.jl",
    format=Documenter.HTML(;
        canonical="https://"Maysam".github.io/ThermalSpreadingResistance.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/"Maysam"/ThermalSpreadingResistance.jl",
    devbranch="main",
)
