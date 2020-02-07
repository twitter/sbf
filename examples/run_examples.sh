cd ..
mvn package
cd examples
java -jar ../target/sbf-1.0.0.jar unweighted8node.config
java -jar ../target/sbf-1.0.0.jar weighted8node.config
