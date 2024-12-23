# For class diagrams
pyreverse -o dot -p peplab peplab/
dot -Grankdir=LR -Tsvg classes_peplab.dot -o classes_peplab.svg

# For package diagrams
pyreverse -o dot -p peplab peplab/
dot -Grankdir=LR -Tsvg packages_peplab.dot -o packages_peplab.svg
