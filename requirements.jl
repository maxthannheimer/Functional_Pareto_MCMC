# Import necessary packages
import Pkg; Pkg.add("Distributions")
import Pkg; Pkg.add("Plots")
import Pkg; Pkg.add("FFTW")
import Pkg; Pkg.add("JLD2")
import Pkg; Pkg.add("LinearAlgebra")
import Pkg; Pkg.add("Printf")
import Pkg; Pkg.add("Random")


# Load required packages
using LinearAlgebra
using Random, Distributions
using Plots
using FFTW
using JLD2
using Printf