
This folder contains python code to simulate the conditioning task in the paper:

Dopamine Pulses Signal Prediction Errors Using A Spectrum Of Time Horizons.  Ali Mohebi1, Wei Wei, Lily Pelattini, Kyoungjun Kim, Joshua D Berke. Nature Neuroscience, 2024.

-----
Packages:

tensorflow. The code was written using tf1.x, but made compatible with tf2 by using tf.compat.v1. 

python 3.8 or higher.

To run:

First run PPO_Meta_ITI.py with a chosen seed or the default one, then run PPO_Meta_ITI_eva.py using a saved checkpoint of the model to do the evaluation.

For questions about the code, please e-mail W. Wei (wei.wei@ucsf.edu)
