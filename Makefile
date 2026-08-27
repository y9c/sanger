# Improved Makefile for sanger
.PHONY: init dev test clean lock

# Install dependencies (production only)
init:
	uv pip install --no-deps

# Install all dependencies (dev + prod)
dev:
	uv pip install

# Run all tests
test:
	python -m unittest discover -s test

# Remove Python cache and temp files
clean:
	rm -rf __pycache__ */__pycache__ *.pyc *.pyo *.pyd temp/* test/__pycache__ sanger/__pycache__

# Update lock file from pyproject.toml
lock:
	uv pip compile
