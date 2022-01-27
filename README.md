# Parallel algorithm for synthetic image generation with application to tokamak plasma diagnostics

The tomographic reconstruction is a powerful diagnostic tool in nuclear fusion experiments for the determination of the shape and position of the plasma. However, neural networks are emerging as a suitable alternative to conventional plasma tomography algorithms. In order to train such AI/ML based models, we need large-scale, diversified image data for learning and evaluation which is difficult to obtain from real experiments. Two different code parallelization strategies have been provided along with the serial algorithm for the generation of this large-scale dataset. The proposed parallel algorithm-II can be used to obtain a large-scale dataset of the synthetic images of plasma within a few seconds which is very important for real time applications.  

# Steps to visualise the synthetic plasma image:

1. Execute any cpp file from the codestack (i.e either Serial.cpp, Algorithm_1.cpp or Algorithm_2.cpp).
2. The code file would create a 'pixels_intensity_txt.txt' file, which would store the intensity value of each pixel generating the synthetic image.
3. Run the 'visualisation.py' file, which reads the 'pixels_intensity_txt.txt' file, and generates the synthetic image of plasma.
4. Input paramaters to the code could be changed in the 'input.txt' file.
