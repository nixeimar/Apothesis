# Apothesis input format

Simulations use JSON version 1. JSON parsing is isolated in
`src/config/configuration.cpp`; the
simulation sees format-independent `Config::Value` and `Config::Configuration`
objects, not JSON library types.

For fixed-site 3D crystal geometry, use `lattice.representation: "atomistic"`
and the [atomistic surface guide](atomistic-surfaces.md). The sections below
describe the existing height-field simulation input.

## Select and validate an input

```sh
./apothesis --input examples/height-field.json
./apothesis -i examples/height-field.json
./apothesis examples/height-field.json
./apothesis --validate-input examples/height-field.json
```

With no arguments, Apothesis opens `input.json` in the working directory.
Every selected input is parsed as JSON. Validation checks the entire JSON input
and referenced matrix files without running the simulation or writing output.
Errors identify the offending field and return a nonzero exit status.
Output files are written in the working directory. Paths to lattice matrix files
are resolved relative to the JSON input file.

Start with the complete [height-field example](../examples/height-field.json). The structural
[JSON Schema](input.schema.json) can be used with an editor or schema validator.
The runtime additionally validates cross references and model restrictions.
Standard JSON does not support comments; use `description` fields for notes.
Duplicate JSON object keys are rejected.

## Structure and extension parameters

| Section | Purpose |
| --- | --- |
| `schema_version` | Required integer `1`; unsupported versions are rejected. |
| `system.parameters` | Temperature, pressure and additional system properties. |
| `simulation` | Timing, random seed and optional simulation parameters. |
| `lattice` | Lattice type and its parameters. |
| `species` | Species definitions keyed by unique names. |
| `processes` | Array of explicitly typed processes with unique IDs. |
| `output` | Output intervals, coverage labels and optional output parameters. |
| `parameters` | Optional root-level application/extension parameters. |

Each species and process, the system, lattice, simulation, output and each rate
model can have a `parameters` object. These objects preserve arbitrary names and
values: null, booleans, signed 64-bit integers, floating-point numbers, strings,
arrays and nested objects. Recognized fields are validated; other parameter values
are retained. Structural misspellings outside parameter objects are rejected.
Custom parameter names inside a parameter object are intentionally unrestricted.

For example, adding `"electric_field": {"strength": 2.0, "axis": [0, 0, 1]}`
to `system.parameters` requires no parser change. It affects the simulation only
when a model reads and uses it. Passing data does not create new physics or register
new process/lattice implementations.

Optional parameter objects are normalized to empty objects in the internal
document so model code can access them consistently.

## System and simulation

`system.parameters.temperature` is positive and measured in K.
`system.parameters.pressure` is nonnegative and measured in Pa.

`simulation.duration` is required, positive and measured in seconds.
`simulation.start_time` defaults to zero. End time is **start_time + duration**.
`simulation.random_seed` is an integer from 0 to 2147483647; absent or zero means
initialize from the current time, preserving existing behavior.

## Lattice and species

`lattice.type` supports `SimpleCubic`, `FCC`, `HCP` and `Diamond`.
Its `parameters` include:

- `dimensions`: `[x, y]`, integer dimensions of at least three sites each. HCP
  requires even x; FCC requires even x and y.
- `initial_species`: a declared species name, used as the unoccupied substrate.
- `initial_height`: a nonnegative integer; or `initial_height_file`, exclusively.
- `initial_species_file`: optional matrix of initial surface labels.
- `growth_species`: optional array of declared species names, without `*`.

Matrix files are supported only for `SimpleCubic`. They contain exactly `y` rows
of `x` whitespace-separated values, without a header. Height values are nonnegative
integers. Species values are declared names, optionally ending in `*` for occupied
surface sites. Missing files, invalid dimensions and undeclared species are errors.
Stepped surfaces are not supported by this JSON version; `steps` is rejected.

Species definitions contain an optional `description` string and `parameters` object.
Use base names such as `CO` or `CO_ads` without a trailing `*`. Names cannot begin
with a digit or dot or contain whitespace or chemical-equation delimiters.
Gas and surface phases are specified on process terms. The current engine receives
`CO` for a gas term and `CO*` for a surface term.

## Processes

A process has a unique `id`, explicit `type`, `reactants`, `products` and `rate`.
An optional `description`, custom `parameters` and boolean `all_neighbors` may
also be supplied. IDs preserve separate channels even when their chemical equations
are identical. The human-readable equation is generated from structured terms;
it is not parsed from an input string or used to infer the process type.

```json
{
  "id": "co_adsorption",
  "type": "adsorption",
  "reactants": [
    {"species": "CO", "phase": "gas"},
    {"vacancy": true}
  ],
  "products": [{"species": "CO", "phase": "surface"}],
  "rate": {"model": "constant", "parameters": {"value": 0.4}},
  "parameters": {"custom_property": [1, 2, 3]}
}
```

Every species reference must be declared. A term has `species`, `phase` (`gas` or
`surface`) and optional positive integer `coefficient` (default 1). A vacancy term
uses `"vacancy": true` and optional `coefficient`, without species or phase.
Combine repeated terms through their coefficient.

Version 1 deliberately validates the existing engine's process shapes:

- `adsorption`: one gas molecule plus one or two vacancies, producing the same
  number of surface species of one kind.
- `desorption`: one surface species producing one vacancy and one gas species.
- `diffusion`: one surface species on each side, with identical names.
- `reaction`: two distinct surface reactant terms, with coefficient one, and one
  or two product terms with coefficient one.

Arbitrary reaction stoichiometries and new physical models require corresponding
engine implementations. `all_neighbors` retains the old `all` behavior for
adsorption, desorption and diffusion. It is unsupported for reactions.

### Rate parameters and units

| Model | Process types | Required named parameters |
| --- | --- | --- |
| `constant` | All four | `value`: nonnegative rate in the existing model's convention; adsorption input uses ML/s, other event rates use s^-1. Existing neighbor factors still apply. |
| `simple` | Adsorption | `sticking_coefficient` and `molar_fraction` in [0, 1]; positive `site_density` in sites/m²; positive `molecular_weight` in kg/mol. |
| `arrhenius` | Reaction, desorption | `prefactor` in s^-1 and `activation_energy` in J/mol, both nonnegative. |
| `arrhenius` | Diffusion | The preceding two parameters plus nonnegative `migration_energy` in J/mol, using the existing neighbor-dependent diffusion expression. |

For `simple`, `molecular_weight` may be omitted from rate parameters when it is
provided in the gas species' `parameters`. A rate-level value takes precedence.
Arrhenius adsorption is not implemented and is rejected.

## Output

`output.log_interval` and `output.lattice_interval` are required, positive values
in seconds. Optional `coverage` uses current engine labels such as `["CO*", "O*",
"X"]`; each base species must be declared. Filenames and output behavior otherwise
follow the existing engine.

## Access custom parameters from C++

```cpp
// From Parameters or a model's getParameters():
const auto& document = parameters->document();
double field = document.at("system").at("parameters")
                       .at("electric_field").at("strength").number();
const auto& species = document.at("species").at("CO").at("parameters");

// Inside a Process subclass; the stable ID survives neighbor expansion:
const auto& custom = inputParameters();
if (custom.contains("enabled")) {
    bool enabled = custom.at("enabled").as<bool>();
    // Use enabled in this model's behavior.
}
```

`number()` accepts integers or floating-point values. Other values use `as<T>()`
with `bool`, `std::int64_t`, `double`, `std::string`, `Config::Value::Array` or
`Config::Value::Object`. Missing keys or incompatible types throw exceptions.
A model must validate custom fields before using them. The full process and rate
objects also remain available through `document().at("processes")`.

The current models still receive positional argument vectors through a small
adapter after named-parameter validation. This preserves their existing interfaces
while allowing models to migrate incrementally to typed parameters.
The JSON library does not leak into model interfaces.

## Tests

CMake builds the parser tests and registers them with CTest. When Python 3 is
available, it also registers executable integration tests:

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

Integration tests cover validation without output, default and explicit JSON selection, explicit process
identity, Arrhenius initialization, and rectangular matrix imports. They run in
temporary directories so simulation outputs do not modify the checkout.
