
.SECONDARY:
.PHONY: clean

CXXFLAGS = -std=c++11

%.o: %.cpp
	$(CXX) -c $< $(CXXFLAGS)

%: %.o
	$(CXX) -o $@ $^ $(CXXFLAGS)

clean:
	rm -rf *.o
