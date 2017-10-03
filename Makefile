SRCS=util.cpp random.cpp block_file.cpp b_node.cpp b_tree.cpp qalsh.cpp ann.cpp main.cpp
OBJS=$(SRCS:.cpp=.o)

CXX?=g++ -std=c++11
CPPFLAGS=-w -O3

.PHONY: clean

all: ${OBJS}
	${CXX} ${CPPFLAGS} -o qalsh ${OBJS}

util.o: util.h

random.o: random.h

block_file.o: block_file.h

b_node.o: b_node.h

b_tree.o: b_tree.h

qalsh.o: qalsh.h

ann.o: ann.h

main.o:

clean:
	-rm ${OBJS} qalsh
