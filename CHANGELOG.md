# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this

## [0.1.5] - 2026-11-26
### Added
- Added OR-Tools support for maximum flow and integer MCCP computations.
- Added MCCP oracle call counts for the QFP.
- Extended the test suite.

### Changed
- Refined the MDFP and QFP implementations to align with the manuscript.
- Improved code readability, documentation, and examples.

### Fixed
- Resolved native-library conflicts between OR-Tools and HiGHS process isolation.
- Improved DIMACS parsing and input validation.

## [0.1.4] - 2026-11-22
### Changed
- Updated the LP based approach to align with the manuscript.

## [0.1.3] - 2026-11-22
### Changed
- Updated the pure interpolation search and safeguarded interpolation search to align with the manuscript.

## [0.1.2] - 2026-11-22
### Changed
- Updated selected methods to use integer precision.

## [0.1.1] - 2026-11-22
### Changed
- Updated the Newton search method to align with the manuscript.

### Added
- Added a validation check to ensure that `F` is an integer in QFP.

## [0.1.0] - 2026-11-05
### Added
- Initialized project.
