# For class diagrams
pyreverse -o dot -p peplab peplab/
dot -Grankdir=LR -Tpdf classes_peplab.dot -o classes_peplab.pdf

# For package diagrams
pyreverse -o dot -p peplab peplab/
dot -Grankdir=LR -Tpdf packages_peplab.dot -o packages_peplab.pdf
