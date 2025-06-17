bld:
	docker build -t crown-gen-server-cpp:1.0 .


run: 
	docker run -p 7812:7812 crown-gen-server-cpp:1.0 