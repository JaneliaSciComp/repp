.PHONY: docs
docs:
	go run ./main.go
	cp ../README.md ./index.md
	bundle install
	bundle exec just-the-docs rake search:init

serve:
	bundle install
	bundle exec jekyll serve
