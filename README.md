# PhyBS

This project is for my PhD independent study, which implements basic and important soft body deformation papers in literature, including force-based simulation, position-based and projective dynamics. Various improvements in numerical solvers and parallelism make it possible to interact with the models in real time. The latest demo shows that for an armadillo model with 43K tets and 12K nodes, the FPS is around 65 running on GPU which is significantly faster than CPU version codes.

The codes was tested with the following configurations:
-Windows 10, Visual Studio 2013, Qt 5.6.2, Eigen 3, CUDA 8.0, Intel Parallel Studio XE 2017.

## Demo
* GPU accelearted elastic body simulation [7,8]
[![Video](https://img.youtube.com/vi/urxVQI-3MAw/0.jpg)](https://www.youtube.com/watch?v=urxVQI-3MAw)

## Rome wasn't built in a day
* FEM/FVM force-based simulation[1]                                                    

[![Video](https://img.youtube.com/vi/RZTT9vSTd5M/0.jpg)](https://www.youtube.com/watch?v=RZTT9vSTd5M)

* Stable simulation using implicit integration[2], accelerated by OpenMP

[![Video](https://img.youtube.com/vi/SMo9IlWolVs/0.jpg)](https://www.youtube.com/watch?v=SMo9IlWolVs)

* Compute stiffness matrix of hyperelasic materials with two different methods: invariants-based[3] and stretch-based[4]

[![Video](https://img.youtube.com/vi/xZR5uczls28/0.jpg)](https://www.youtube.com/watch?v=xZR5uczls28)

* Position based dynamics using strain constraints[5]                              

[![Video](https://img.youtube.com/vi/HgDR9nFfIRs/0.jpg)](https://www.youtube.com/watch?v=HgDR9nFfIRs)


## References
[1] Teran, J., Blemker, S., Ng Thow Hing, V. and Fedkiw, R., Finite Volume Methods for the Simulation of Skeletal Muscle, ACM SIGGRAPH/Eurographics Symposium on Computer Animation (SCA), edited by D. Breen and M. Lin, pp. 68-74, 2003.

[2] David B. , Andrew W. 1998. Large steps in cloth simulation. In Proceedings of the 25th annual conference on Computer graphics and interactive techniques (SIGGRAPH '98). ACM, New York, NY, USA, 43-54

[3] Teran, J., Sifakis, E., Irving, G. and Fedkiw, R., "Robust Quasistatic Finite Elements and Flesh Simulation", ACM SIGGRAPH/Eurographics Symposium on Computer Animation (SCA), edited by K. Anjyo and P. Faloutsos, pp. 181-190, 2005.

[4] Hongyi X. , Fun Shing S. , Yufeng Z. , Jernej B. Nonlinear Material Design Using Principal Stretches, ACM Transactions on Graphics 34(4) (SIGGRAPH 2015), Los Angeles, CA, USA

[5] Müller, M., Chentanez, N., Kim, T. Y., & Macklin, M. (2014, July). Strain based dynamics. In Proceedings of the ACM SIGGRAPH/Eurographics Symposium on Computer Animation (pp. 149-157). Eurographics Association.

[6] Bouaziz, S., Martin, S., Liu, T., Kavan, L., & Pauly, M. (2014). Projective dynamics: fusing constraint projections for fast simulation. ACM Transactions on Graphics (TOG), 33(4), 154.

[7] Wang, H. (2015). A chebyshev semi-iterative approach for accelerating projective and position-based dynamics. ACM Transactions on Graphics (TOG), 34(6), 246.

[8] Wang, H., & Yang, Y. (2016). Descent methods for elastic body simulation on the GPU. ACM Transactions on Graphics (TOG), 35(6), 212.