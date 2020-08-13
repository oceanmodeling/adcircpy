default: error

coverage:
	@coverage run --source adcircpy -m nose --rednose --nologcapture --verbose && coverage report -m && coverage-badge -f -o tests/coverage.svg && rm -rf .coverage

badge:
	@$(shell git diff --exit-code -s tests/coverage.svg)
	@if ! [ $(.SHELLSTATUS) = 0 ]; then git add tests/coverage.svg && git commit -m "Updated coverage badge."; fi

check-tag:
ifndef tag
	$(error tag is undefined)
endif

pre-release: check-tag coverage
	@git tag $(tag); rm -rf ./dist; ./setup.py sdist;
	
release: pre-release
	@git push; git push --tags; twine upload dist/*

error:
	@echo -e "ERROR: Argument required."
	@exit 2