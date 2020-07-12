default: error

coverage:
	@coverage run --source adcircpy -m nose --rednose --nologcapture --verbose && coverage report -m && coverage-badge -f -o tests/coverage.svg && rm -rf .coverage

error:
	@echo -e "ERROR: Argument required."
	@exit 2
