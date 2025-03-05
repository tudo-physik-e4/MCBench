"""
    abstract type AnySampler

An abstract type that serves as a base for all sampling algorithms.
"""
abstract type AnySampler end
"""
    abstract type SamplingAlgorithm <: AnySampler

An abstract type for general sampling algorithms that inherit from `AnySampler`.
"""
abstract type SamplingAlgorithm  <: AnySampler end
"""
    abstract type IIDSamplingAlgorithm <: SamplingAlgorithm

An abstract type for independent and identically distributed (IID) sampling algorithms that inherit from `SamplingAlgorithm`.
This type doesn't have any fields or methods, but it is used to have testfunctions reference their IID `Base.rand` method when sampling. 
"""
abstract type IIDSamplingAlgorithm <: SamplingAlgorithm end
"""
    abstract type AbstractFileBasedSampler <: AnySampler

An abstract type for file-based sampling algorithms that inherit from `AnySampler`.
"""
abstract type AbstractFileBasedSampler <: AnySampler end

"""
    struct IIDSampler <: IIDSamplingAlgorithm

A struct representing an IID (Independent and Identically Distributed) Sampler to created instances of IID sampling algorithms.

# Fields
- `n_steps::Int`: The number of steps for the sampler. Identical to the number of samples for IID.
- `info::String`: Information or description of the sampler.

# Constructors
- `IIDSampler()`: Creates an `IIDSampler` with default values of `10^5` steps and "IID" as the info string.
"""
struct IIDSampler <: IIDSamplingAlgorithm 
    n_steps::Int
    info::String
end

IIDSampler() = IIDSampler(10^5,"IID")

"""
    struct FileBasesSampler <: SamplingAlgorithm
IO funtionalities to read data from a set of files.

# Fields
- `files::Vector{String}`: Vector of file paths.
- `current_file_index::Int`: Index of the current file in the vector.
- `current_position::Int`: Position within the current file.
- `current_file_handle::IOStream`: Handle to the currently open file.
- `info::String`: Information about the sampler. Used for plotting.

# Constructors
- `FileBasedSampler(; fields...)`
- `FileBasedSampler(file_paths::Vector{String})`: Creates a `FileBasedSampler` with the given vector of file paths.
- `FileBasedSampler(path::String)`: Creates a `FileBasedSampler` with the given file path. This will load all files in the directory if the path is a directory.
"""
mutable struct FileBasedSampler <: AbstractFileBasedSampler
    files::Vector{String}           # Vector of file paths
    current_file_index::Int         # Index of the current file in the vector
    current_position::Int           # Position within the current file
    current_file_handle::IOStream   # Handle to the currently open file
    info::String                    # Information about the sampler
end

function FileBasedSampler(file_paths::Vector{String}; info::String="FileBasedSampler")
    @assert !isempty(file_paths) "File paths cannot be empty"
    # Open the first file
    first_file_handle = open(file_paths[1], "r")
    return FileBasedSampler(file_paths, 1, 1, first_file_handle, info)
end

function FileBasedSampler(path::String; info::String="FileBasedSampler")
    # Determine if the path is a file or a directory
    if isfile(path)
        files = [path]  # Single file
    elseif isdir(path)
        files = sort(readdir(path, join=true))  # All files in the directory
    else
        error("Path does not exist or is neither a file nor a directory: $path")
    end

    @assert !isempty(files) "No files to read from."
    first_file_handle = open(files[1], "r")
    return FileBasedSampler(files, 1, 1, first_file_handle, info)
end

function read_sample!(sampler::FileBasedSampler)
    # If the current file handle is at the end, move to the next file
    if eof(sampler.current_file_handle)
        close(sampler.current_file_handle)
        sampler.current_file_index += 1
        if sampler.current_file_index > length(sampler.files)
            error("No more files to read from.")
        end
        sampler.current_file_handle = open(sampler.files[sampler.current_file_index], "r")
        sampler.current_position = 0
    end

    # Read the next line/sample
    sample = readline(sampler.current_file_handle)
    sampler.current_position += 1
    return sample
end

function reset_sampler!(sampler::FileBasedSampler)
    # Reset to the first file and position
    close(sampler.current_file_handle)
    sampler.current_file_index = 1
    sampler.current_position = 1
    sampler.current_file_handle = open(sampler.files[1], "r")
end

function close_sampler!(sampler::FileBasedSampler)
    close(sampler.current_file_handle)
end

"""
    struct CsvBasedSampler <: AbstractFileBasedSampler
    IO funtionalities to read data from a set of CSV files.

    # Fields
    - `fbs::FileBasedSampler`: File-based sampler
    - `header::Vector{String}`: Header of the CSV file
    - `mask::Vector{Int}`: Mask to extract the desired columns
    - `info::String`: Information about the sampler

    # Constructors
    - `CsvBasedSampler(; fields...)`
    - `CsvBasedSampler(file_paths::Vector{String})`: Creates a `CsvBasedSampler` with the given vector of file paths.
    - `CsvBasedSampler(path::String)`: Creates a `CsvBasedSampler` with the given file path. This will load all files in the directory if the path is a directory.
"""
mutable struct CsvBasedSampler <: AbstractFileBasedSampler
    fbs::FileBasedSampler           # File-based sampler
    header::Vector{String}          # Header of the CSV file
    mask::Vector{Int}               # Mask to extract the desired columns
    info::String                    # Information about the sampler
end

function CsvBasedSampler(file_paths::Vector{String}; info::String="CsvBasedSampler")
    fbs = FileBasedSampler(file_paths,info=info)
    header = split(read_sample!(fbs), ",")
    return CsvBasedSampler(fbs,header, collect(1:length(header)),info)
end

function CsvBasedSampler(path::String; info::String="CsvBasedSampler")
    fbs = FileBasedSampler(path,info=info)
    header = split(read_sample!(fbs), ",")
    return CsvBasedSampler(fbs,header, collect(1:length(header)),info)
end

function read_sample!(sampler::CsvBasedSampler)
    i = sampler.fbs.current_file_handle
    sm = read_sample!(sampler.fbs) #split(read_sample!(sampler.fbs),",")
    j = sampler.fbs.current_file_handle
    if i == j       #if next file skip header line
        # return [parse(Float64,i) for i in sm[sampler.mask]]
        return sm
    else
        # return [parse(Float64,i) for i in split(read_sample!(sampler.fbs),",")[sampler.mask]]
        return read_sample!(sampler.fbs)
    end
end

function set_mask(sampler::CsvBasedSampler, seq::Vector{String})
    seq_symbols = Symbol.(seq) 
    header_map = Dict(Symbol(h) => i for (i, h) in enumerate(sampler.header))
    sampler.mask = [header_map[col] for col in seq_symbols]
end

function reset_sampler!(sampler::CsvBasedSampler)
    reset_sampler!(sampler.fbs)
    header = split(read_sample!(sampler.fbs), ",")
end

""" 
    struct DsvSampler{D<:DensitySampleVector} <: AbstractFileBasedSampler
    IO funtionalities to read data from a set of DensitySampleVector files.

    # Fields
    - `dsvs::Vector{D}`: Vector of file paths
    - `current_dsv_index::Int`: Index of the current file in the vector
    - `current_position::Int`: Position within the current file
    - `weighted::Bool`: Whether the samples are weighted
    - `total_samples::Int`: Total number of samples
    - `neff::Vector{Float64}`: Number of effective samples
    - `info::String`: Information about the sampler

    # Constructors
    - `DsvSampler(; fields...)`
    - `DsvSampler(dsvs::Vector{D})`: Creates a `DsvSampler` with the given vector of file paths.
"""
mutable struct DsvSampler{D<:DensitySampleVector} <: AbstractFileBasedSampler
    dsvs::Vector{D}   # Vector of file paths
    current_dsv_index::Int              # Index of the current file in the vector
    current_position::Int               # Position within the current file
    weighted::Bool                      # Whether the samples are weighted
    total_samples::Int                  # Total number of samples
    neff::Vector{Float64}               # Number of effective samples
    info::String                        # Information about the sampler
end

function DsvSampler(dsvs::Vector{D}; info::String="DsvSampler") where {D<:DensitySampleVector}
    @assert !isempty(dsvs) "Density sample vectors cannot be empty"
    weighted = false
    total_samples = 0
    neff = Float64[]
    for dsv in dsvs
        total_samples += sum(dsv.weight)
        push!(neff, get_effective_sample_size(dsv))
        if dsv.weight != ones(length(dsv))
            weighted = true
            #break
        end
    end
    return DsvSampler(dsvs, 1, 1, weighted, Int(total_samples), neff, info)
end

function read_sample!(sampler::DsvSampler)
    if sampler.current_position > length(sampler.dsvs[sampler.current_dsv_index])
        sampler.current_dsv_index += 1
        sampler.current_position = 1
    end
    # If the current file handle is at the end, move to the next file
    if sampler.current_dsv_index > length(sampler.dsvs)
        error("No more files to read from.")
    end
    # Read the next line/sample
    sample = sampler.dsvs[sampler.current_dsv_index][sampler.current_position]
    sampler.current_position += 1
    return sample
end

function unweight!(sampler::DsvSampler)
    # for i in 1:length(sampler.dsvs)
    #     dsv_unw = sampler.dsvs[i]
    #     dsv_unw = resample_dsv_to_ess(dsv_unw)
    #     sampler.dsvs[i] = dsv_unw
    # end
end