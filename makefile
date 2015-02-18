
.SECONDARY:
.PHONY: clean

%.o: %.cpp
	$(CXX) -c $< $(CXXFLAGS)

clean:
	rm -rf *.o
