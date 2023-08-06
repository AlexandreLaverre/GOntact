build:
	dune build
	dune build server/gontact_client.bc.js
	mkdir -p server/static/js
	rm -f server/static/js/gontact_client.js
	cp _build/default/server/gontact_client.bc.js server/static/js/gontact_client.js

%.exe:
	dune build $@

local: build
	dune exec -- gontact_server
