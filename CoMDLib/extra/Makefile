all: libDriver libVerifier

libDriver: libDriver.c libCoMD.a
	cc -lm -lCoMD -L./  -lCoMD libDriver.c -lCoMD -o ./libDriver

libVerifier: libVerifier.c libCoMD.a
	cc -lm -lCoMD -L./  -lCoMD libVerifier.c -lCoMD -o ./libVerifier

clean:
	rm -f ./libDriver ./libVerifier

