pyreverse -o dot -p peplab peplab/
dot -Grankdir=LR -Tpng classes_peplab.dot -o classes_peplab.png
pyreverse -o dot -p peplab peplab/
dot -Grankdir=LR -Tpng packages_peplab.dot -o packages_peplab.png
