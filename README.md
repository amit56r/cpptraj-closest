 
 
 The following need to be done in order to integrate into cpptraj 
 
 - makefile should use nvcc -arch=sm_xx instead of gcc/g++
 - The following code/def needs to be added to action_closest.cpp/action_closest.h
 ```C++
    void Action_NoImage(Frame& frmIn,double maxD);
    void Action_ImageOrtho(Frame& frmIn, double maxD);
    void Action_ImageNonOrtho(Frame& frmIn, double maxD, Matrix_3x3 ucell, Matrix_3x3 recip);

    bool cuda_action_center(Frame& frmIn, double maxD, Matrix_3x3 ucell, Matrix_3x3 recip,int type, float &time_gpu);
    bool cuda_action_no_center(Frame& frmIn, double maxD, Matrix_3x3 ucell, Matrix_3x3 recip,int type, float &time_gpu);
  ```
    
- The folder `Cuda_kernels` needs to be copied
