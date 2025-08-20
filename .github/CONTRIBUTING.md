## Contributing to CRISPR Toolkit

Thank you for your interest in contributing to the CRISPR Toolkit! This document provides guidelines and information about contributing to this project.

### Code of Conduct

This project and everyone participating in it is governed by our Code of Conduct. By participating, you are expected to uphold this code.

### How Can I Contribute?

#### Reporting Bugs

Before creating bug reports, please check the existing issues as you might find out that you don't need to create one. When you are creating a bug report, please include as many details as possible:

- Use a clear and descriptive title
- Describe the exact steps which reproduce the problem
- Provide specific examples to demonstrate the steps
- Describe the behavior you observed after following the steps
- Explain which behavior you expected to see instead and why
- Include details about your configuration and environment

#### Suggesting Enhancements

Enhancement suggestions are tracked as GitHub issues. When creating an enhancement suggestion, please include:

- Use a clear and descriptive title
- Provide a step-by-step description of the suggested enhancement
- Provide specific examples to demonstrate the steps
- Describe the current behavior and explain which behavior you expected to see instead
- Explain why this enhancement would be useful

#### Pull Requests

Please follow these steps to have your contribution considered by the maintainers:

1. Follow all instructions in the template
2. Follow the styleguides
3. After you submit your pull request, verify that all status checks are passing

### Development Process

1. Fork the repository
2. Create a new branch from `develop` for your feature or bug fix
3. Make your changes
4. Add or update tests as necessary
5. Ensure all tests pass
6. Update documentation as necessary
7. Submit a pull request

### Style Guidelines

#### Python Style Guide

- Follow PEP 8
- Use Black for code formatting
- Use type hints where appropriate
- Write docstrings for all functions, classes, and modules
- Maximum line length: 88 characters

#### Java Style Guide

- Follow Google Java Style Guide
- Use meaningful variable and method names
- Write Javadoc for public methods and classes
- Use camelCase for variable and method names
- Use PascalCase for class names

#### C++ Style Guide

- Follow Google C++ Style Guide
- Use snake_case for variable and function names
- Use PascalCase for class names
- Use ALL_CAPS for constants
- Include proper header guards

#### Git Commit Messages

- Use the present tense ("Add feature" not "Added feature")
- Use the imperative mood ("Move cursor to..." not "Moves cursor to...")
- Limit the first line to 72 characters or less
- Reference issues and pull requests liberally after the first line

### Testing

- Write tests for all new functionality
- Ensure existing tests continue to pass
- Aim for high test coverage
- Use pytest for Python tests
- Use JUnit for Java tests
- Use Google Test for C++ tests

### Documentation

- Update README.md if necessary
- Update relevant documentation in the `docs/` directory
- Include docstrings/comments for new code
- Update CHANGELOG.md for significant changes

### License

By contributing to CRISPR Toolkit, you agree that your contributions will be licensed under the MIT License.
