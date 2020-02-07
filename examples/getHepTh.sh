wget http://konect.uni-koblenz.de/downloads/tsv/ca-cit-HepTh.tar.bz2
bunzip2 ca-cit-HepTh.tar.bz2
tar -zxvf ca-cit-HepTh.tar
cut -f1,2 -d" " ca-cit-HepTh/out.ca-cit-HepTh | sed '/^%/d' | sort | uniq > hepTh.edges


