#include "BSplinePose.hpp"
#include "assert_macros.hpp"
// boost::tie
#include <boost/tuple/tuple.hpp>
#include "transformations.hpp"

namespace bsplines {
  using namespace sm::kinematics;

    // given T, T'\in\mathrm{SE}(3), define Exp_{T}(h), h\in\mathbb{R}^6 , Log_{T}(T') 
    // [p,\theta]\in\mathbb{R}^6 , \mathrm{SE}(3) \ni T = Exp[p, \theta] = [exp(-\theta) p]
    //                                                                     [      0      1] 
    // T = [R t] \in\mathrm{SE}(3), Log T = [t  -log(R)]
    //     [0 1]
    // T \oplus [a, b] = Exp_{T}([a, b]) = Exp[a, b] * T
    // T1 \ominus T2 = Log_{T2}(T1) = Log(T1 * T2^{-1})
    BSplinePose::BSplinePose(int splineOrder, const RotationalKinematics::Ptr & rotationalKinematics) 
      : BSpline(splineOrder), rotation_(rotationalKinematics)
    {
      
    }

    BSplinePose::~BSplinePose()
    {

    }
      
    Eigen::Matrix4d BSplinePose::transformation(double tk) const
    {
      return curveValueToTransformation(eval(tk));
    }

    // T: transformation  p: axis angle  C: spline control points
    // dT/dC = dT/dp * dp/dC
    Eigen::Matrix4d BSplinePose::transformationAndJacobian(double tk, Eigen::MatrixXd * J, Eigen::VectorXi * coefficientIndices) const
    {
      Eigen::MatrixXd JS;  // dp/dC
      Eigen::VectorXd p;
      p = evalDAndJacobian(tk,0,&JS, coefficientIndices);
      
      Eigen::MatrixXd JT;  // dT/dp
      Eigen::Matrix4d T = curveValueToTransformationAndJacobian( p, &JT );      
      
      if(J)
      {
        *J = JT * JS;
      }

      return T;
    }

    // Log_{O(\theta)}O(\theta + \delta) = J\delta
    // => J\delta = Log(O(\theta + \delta) * O(\theta)^{-1})
    //            = Log(exp(-(\theta + \delta)) * exp(\theta))
    //            = Log(exp(-Jr(\theta)\delta))
    //            = Jr(\theta)\delta
    Eigen::Matrix3d BSplinePose::orientationAndJacobian(double tk, Eigen::MatrixXd * J, Eigen::VectorXi * coefficientIndices) const
    {
      Eigen::Matrix3d C;
      Eigen::MatrixXd JS;
      Eigen::VectorXd p;
      p = evalDAndJacobian(tk,0,&JS, coefficientIndices);

      Eigen::Matrix3d S;
      C = rotation_->parametersToRotationMatrix(p.tail<3>(), &S);

      Eigen::MatrixXd JO = Eigen::MatrixXd::Zero(3,6);
      JO.block(0,3,3,3) = S;
      if(J)
      {
        *J = JO * JS;
      }

      return C;

    }

    // Log_{IO(\theta)}IO(\theta + \delta) = J\delta
    // => J\delta = Log(IO(\theta + \delta) * IO(\theta)^{-1})
    //            = Log(exp(\theta + \delta) * exp(-\theta))
    //            = Log(Jl(\theta)\delta)
    //            = -Jl(\theta)\delta
    //            = -exp(\theta)Jr(\theta)\delta
    Eigen::Matrix3d BSplinePose::inverseOrientationAndJacobian(double tk, Eigen::MatrixXd * J, Eigen::VectorXi * coefficientIndices) const
    {
      Eigen::Matrix3d C;
      Eigen::MatrixXd JS;
      Eigen::VectorXd p;
      p = evalDAndJacobian(tk,0,&JS, coefficientIndices);

      Eigen::Matrix3d S;
      C = rotation_->parametersToRotationMatrix(p.tail<3>(), &S).transpose();

      Eigen::MatrixXd JO = Eigen::MatrixXd::Zero(3,6);
      JO.block(0,3,3,3) = S;
      if(J)
      {
        *J = -C * JO * JS;
      }

      return C;

    }


    Eigen::Matrix4d BSplinePose::inverseTransformationAndJacobian(double tk, Eigen::MatrixXd * J, Eigen::VectorXi * coefficientIndices) const
    {
      //std::cout << __FUNCTION__ << "()\n";
      //ASRL_THROW(std::runtime_error,"Not Implemented");
      Eigen::MatrixXd JS;
      Eigen::VectorXd p;
      p = evalDAndJacobian(tk,0,&JS, coefficientIndices);
      
      Eigen::MatrixXd JT;
      Eigen::Matrix4d T = curveValueToTransformationAndJacobian( p, &JT );   
      // Invert the transformation.
      T.topLeftCorner<3,3>().transposeInPlace();
      T.topRightCorner<3,1>() = (-T.topLeftCorner<3,3>() * T.topRightCorner<3,1>()).eval();      

      if(J)
      {
        // The "box times" is the linearized transformation way of inverting the jacobian.
        *J = -sm::kinematics::boxTimes(T) * JT * JS;
      }

      if(coefficientIndices)
      {
        *coefficientIndices = localCoefficientVectorIndices(tk);
      }
      
      return T;
    }
    
    Eigen::Matrix4d BSplinePose::inverseTransformation(double tk) const
    {
      Eigen::Matrix4d T = curveValueToTransformation(eval(tk));
      T.topLeftCorner<3,3>().transposeInPlace();
      T.topRightCorner<3,1>() = (-T.topLeftCorner<3,3>() * T.topRightCorner<3,1>()).eval();
      return T;
    }


    // d(Tv)/dC = d(Tv)/dT * dT/dp * dp/dC
    // d(Tv)/dT s.t. Exp_{T}(h)v - Tv = d(Tv)/dT * h
    Eigen::Vector4d BSplinePose::transformVectorAndJacobian(double tk, const Eigen::Vector4d & v_tk, Eigen::MatrixXd * J, Eigen::VectorXi * coefficientIndices) const
    {
      Eigen::MatrixXd JT;
      Eigen::Matrix4d T_n_vk = transformationAndJacobian(tk, &JT, coefficientIndices);
      Eigen::Vector4d v_n = T_n_vk * v_tk;
      
      if(J)
      {
        *J = sm::kinematics::boxMinus(v_n) * JT;
      }

      return v_n;
    }

    Eigen::Vector4d BSplinePose::inverseTransformVectorAndJacobian(double tk, const Eigen::Vector4d & v_tk, Eigen::MatrixXd * J, Eigen::VectorXi * coefficientIndices) const
    {
      Eigen::MatrixXd JT;
      Eigen::Matrix4d T_n_vk = inverseTransformationAndJacobian(tk, &JT, coefficientIndices);
      Eigen::Vector4d v_n = T_n_vk * v_tk;

      if(J)
      {
        *J = sm::kinematics::boxMinus(v_n) * JT;
      }

      return v_n;
    }
      
    Eigen::Vector3d BSplinePose::position(double tk) const
    {
      return eval(tk).head<3>();
    }



    Eigen::Matrix3d BSplinePose::orientation(double tk) const
    {
      return rotation_->parametersToRotationMatrix(eval(tk).tail<3>());
    }

    Eigen::Matrix3d BSplinePose::inverseOrientation(double tk) const
    {
      return rotation_->parametersToRotationMatrix(eval(tk).tail<3>()).transpose();
    }



    Eigen::Vector3d BSplinePose::linearVelocity(double tk) const
    {
      return evalD(tk,1).head<3>();
    }

    Eigen::Vector3d BSplinePose::linearVelocityBodyFrame(double tk) const
    {
      Eigen::VectorXd r = evalD(tk, 0);
      Eigen::Matrix3d C_wb = rotation_->parametersToRotationMatrix(r.tail<3>());
      return C_wb.transpose() * evalD(tk, 1).head<3>();
    }

    Eigen::Vector3d BSplinePose::linearAcceleration(double tk) const
    {
      return evalD(tk,2).head<3>();
    }

    Eigen::Vector3d BSplinePose::linearAccelerationBodyFrame(double tk) const
    {
      Eigen::VectorXd r = evalD(tk, 0);
      Eigen::Matrix3d C_wb = rotation_->parametersToRotationMatrix(r.tail<3>());
      return C_wb.transpose() * evalD(tk, 2).head<3>();
    }

    Eigen::Vector3d BSplinePose::linearAccelerationAndJacobian(double tk, Eigen::MatrixXd * J, Eigen::VectorXi * coefficientIndices) const
    {
      
      Eigen::Vector3d a = evalDAndJacobian(tk,2,J,coefficientIndices).head<3>();
      if(J)
      {
        J->conservativeResize(3,J->cols());
      }
      return a;
    }
    
    // {}_w\omege_{wb}
    // \omega_w_{b,w} (angular velocity of the body frame as seen from the world frame, expressed in the world frame)
    Eigen::Vector3d BSplinePose::angularVelocity(double tk) const
    {
      Eigen::Vector3d omega;
      Eigen::VectorXd r = evalD(tk,0);
      Eigen::VectorXd v = evalD(tk,1);

      omega = -rotation_->parametersToSMatrix(r.tail<3>()) * v.tail<3>();
      return omega;
    }

    // {}_b\omege_{wb}
    // \omega_w_{b,b} (angular velocity of the body frame as seen from the world frame, expressed in the body frame)
    Eigen::Vector3d BSplinePose::angularVelocityBodyFrame(double tk) const
    {
      Eigen::Vector3d omega;
      Eigen::VectorXd r = evalD(tk,0);
      Eigen::VectorXd v = evalD(tk,1);
      Eigen::Matrix3d S;
      Eigen::Matrix3d C_w_b = rotation_->parametersToRotationMatrix(r.tail<3>(), &S);

      omega = -C_w_b.transpose() * S * v.tail<3>();
      return omega;      

    }

    // \omega_w_{b,b} (angular velocity of the body frame as seen from the world frame, expressed in the body frame)
    Eigen::Vector3d BSplinePose::angularVelocityBodyFrameAndJacobian(double tk, Eigen::MatrixXd * J, Eigen::VectorXi * coefficientIndices) const
    {
      Eigen::MatrixXd Jr;
      Eigen::Matrix3d C_b_w = inverseOrientationAndJacobian(tk,&Jr,NULL);
      
      Eigen::Vector3d omega = angularVelocityAndJacobian(tk, J, coefficientIndices);
      omega = C_b_w * omega;

      if(J)
      {
        *J = C_b_w * (*J) + sm::kinematics::crossMx(omega) * Jr;
      }

      return omega;
    }


    // \omega_w_{b,w} (angular velocity of the body frame as seen from the world frame, expressed in the world frame)
    Eigen::Vector3d BSplinePose::angularVelocityAndJacobian(double tk, Eigen::MatrixXd * J, Eigen::VectorXi * coefficientIndices) const
    {

      Eigen::Vector3d omega;
      Eigen::Vector3d p;
      Eigen::Vector3d pdot;
      Eigen::MatrixXd Jp;
      Eigen::MatrixXd Jpdot;
      p = evalDAndJacobian(tk,0,&Jp,NULL).tail<3>();
      pdot = evalDAndJacobian(tk,1,&Jpdot,coefficientIndices).tail<3>();
      
      // Rearrange the spline jacobian matrices. Now Jpdot is the 
      // jacobian of p wrt the spline coefficients stacked on top
      // of the jacobian of pdot wrt the spline coefficients.
      Jpdot.block(0,0,3,Jpdot.cols()) = Jp.block(3,0,3,Jp.cols());
      
      //std::cout << "Jpdot\n" << Jpdot << std::endl;
      // \omega = -Jr(\theta) * \dot{\theta}
      // \theta, \dot{\theta} are independent, then 
      //    d\omega/dC = d\omega/d\theta * d\theta/dC + d\omega/d\dot{\theta} * d\dot{\theta}/dC
      //  d\omega/d\theta * d\theta/dC + d\omega/d\dot{\theta} * d\dot{\theta}/dC
      // = -R * [d(Jr(\theta)\dot{\theta})/d\theta   Jr(\theta)] * [    d\theta/dC    ]
      //                                                           [ d\dot{\theta}/dC ]
      Eigen::Matrix<double,3,6> Jo;
      omega = -rotation_->angularVelocityAndJacobian(p,pdot,&Jo);
      
      //std::cout << "Jo:\n" << Jo << std::endl;
      if(J)
      {
        *J = -Jo * Jpdot;
      }

      return omega;
    }

    // {}_b\dot\omega_{wb}
    // \omega_dot_w_{b,b} (angular acceleration of the body frame as seen from the world frame, expressed in the body frame)
    Eigen::Vector3d BSplinePose::angularAccelerationBodyFrame(double tk) const
    {
    	Eigen::Vector3d p = evalD(tk,0).tail<3>();
      Eigen::Vector3d v = evalD(tk,1).tail<3>();
    	Eigen::Vector3d a = evalD(tk,2).tail<3>();
    	Eigen::Matrix3d C_w_b = rotation_->parametersToRotationMatrix(p);

      Eigen::Matrix<double,3,6> Jo;
      rotation_->angularVelocityAndJacobian(p,v,&Jo);

      Eigen::Matrix<double,6,1> va;
      va << v, a;

      return -C_w_b.transpose() * (Jo * va);
    }

    // \omega_dot_w_{b,b} (angular acceleration of the body frame as seen from the world frame, expressed in the body frame)
    Eigen::Vector3d BSplinePose::angularAccelerationBodyFrameAndJacobian(double tk, Eigen::MatrixXd * J, Eigen::VectorXi * coefficientIndices) const
    {
      Eigen::MatrixXd Jr;
    	Eigen::Matrix3d C_b_w = inverseOrientationAndJacobian(tk,&Jr,NULL);

      Eigen::Vector3d acce = angularAccelerationAndJacobian(tk, J, coefficientIndices);
      acce = C_b_w * acce;

      if(J)
    	{
    		*J = C_b_w * (*J) + sm::kinematics::crossMx(acce) * Jr;
    	}

      return acce;
    }

    // {}_w\dot\omega_{wb}
    // \omega_dot_w_{b,w} (angular acceleration of the body frame as seen from the world frame, expressed in the world frame)
    Eigen::Vector3d BSplinePose::angularAcceleration(double tk) const
    {
    	Eigen::Vector3d p = evalD(tk,0).tail<3>();
      Eigen::Vector3d v = evalD(tk,1).tail<3>();
    	Eigen::Vector3d a = evalD(tk,2).tail<3>();

      Eigen::Matrix<double,3,6> Jo;
      rotation_->angularVelocityAndJacobian(p,v,&Jo);

      Eigen::Matrix<double,6,1> va;
      va << v, a;

      return -(Jo * va);
    }

    // \todo Only support RotationVector by now
    // \omega_dot_w_{b,w} (angular acceleration of the body frame as seen from the world frame, expressed in the world frame)
    Eigen::Vector3d BSplinePose::angularAccelerationAndJacobian(double tk, Eigen::MatrixXd * J, Eigen::VectorXi * coefficientIndices) const
    {
      Eigen::Vector3d p;
      Eigen::Vector3d v;
    	Eigen::Vector3d a;
      Eigen::MatrixXd Jp;
      Eigen::MatrixXd Jv;
      Eigen::MatrixXd Ja;
      p = evalDAndJacobian(tk,0,&Jp,NULL).tail<3>();
      v = evalDAndJacobian(tk,1,&Jv,NULL).tail<3>();
      a = evalDAndJacobian(tk,2,&Ja,coefficientIndices).tail<3>();

      Eigen::Matrix<double,3,6> Jo;
      rotation_->angularVelocityAndJacobian(p,v,&Jo);

      Eigen::Matrix<double,6,1> va;
      va << v, a;

      Eigen::Vector3d acce = -(Jo * va);

      // Rearrange the spline jacobian matrices. Now Ja is the
    	// jacobian of p wrt the spline coefficients stacked on top
    	// of the jacobian of a wrt the spline coefficients.
    	Ja.block(0,0,3,Ja.cols()) = Jp.block(3,0,3,Jp.cols());

      // Rearrange the spline jacobian matrices. Now Jv is the
    	// jacobian of p wrt the spline coefficients stacked on top
    	// of the jacobian of v wrt the spline coefficients.
    	Jv.block(0,0,3,Jv.cols()) = Jp.block(3,0,3,Jp.cols());

      rotation_->angularVelocityAndJacobian(p,a,&Jo);

      // f(\theta, \dot\theta) = \frac{\partial Jr(\theta)\dot\theta}{\partial\theta}\dot\theta
      // [ \frac{\partial f}{\partial\theta}  \frac{\partial f}{\partial\dot\theta} ]
      Eigen::Matrix<double,3,6> Jpp;

      double factor[5];
      const double pv = p.dot(v);
      const double v2 = v.squaredNorm();
      const double p2 = p.squaredNorm();
      if(p2 < 1e-14)
      {
        const double p4 = p2 * p2;
        // Series[(x Sin[x] + 2 Cos[x] - 2) / x^4, {x, 0, 4}]
        factor[0] = -1.0 / 12 + p2 / 180 - p4 / 6720;
        // Series[(x - Sin[x]) / x^3, {x, 0, 4}]
        factor[1] = 1.0 / 6 - p2 / 120 + p4 / 5040;
        // Series[(3 Sin[x] - x Cos[x] - 2 x) / (x^5), {x, 0, 4}]
        factor[2] = -1.0 / 60 + p2 / 1260 - p4 / 60480;
        // Series[((x^2 - 8) Cos[x] - 5 x Sin[x] + 8) / x^6, {x, 0, 4}]
        factor[3] = 1.0 / 90 - p2 / 1680 + p4 / 75600;
        // Series[((x^2 - 15) Sin[x] + 7 x Cos[x] + 8 x) / x^7, {x, 0, 4}]
        factor[4] = 1.0 / 630 - p2 / 15120 + p4 / 831600;
      }
      else
      {
        const double p1 = std::sqrt(p2);
        const double p3 = p1 * p2;
        const double p4 = p1 * p3;
        const double p5 = p1 * p4;
        const double p6 = p1 * p5;
        const double p7 = p1 * p6;
        const double cosp = std::cos(p1);
        const double sinp = std::sin(p1);

        factor[0] = (p1 * sinp + 2 * cosp - 2) / p4;
        factor[1] = (p1 - sinp) / p3;
        factor[2] = (3 * sinp - p1 * cosp - 2 * p1) / p5;
        factor[3] = ((p2 - 8) * cosp - 5 * p1 * sinp + 8) / p6;
        factor[4] = ((p2 - 15) * sinp + 7 * p1 * cosp + 8 * p1) / p7;
      }

      Jpp.leftCols<3>() = factor[0] * sm::kinematics::crossMx(v) * (pv * Eigen::Matrix3d::Identity() + p * v.transpose()) +
                          factor[1] * (v2 * Eigen::Matrix3d::Identity() - v * v.transpose()) +
                          factor[2] * (pv * pv * Eigen::Matrix3d::Identity() + 2 * pv * p * v.transpose() - 2 * pv * v * p.transpose() - p2 * v * v.transpose()) +
                          factor[3] * pv * sm::kinematics::crossMx(v) * p * p.transpose() +
                          factor[2] * (v2 * p * p.transpose() - pv * v * p.transpose()) +
                          factor[4] * (pv * pv * p * p.transpose() - pv * p2 * v * p.transpose());
      Jpp.rightCols<3>() = factor[0] * (sm::kinematics::crossMx(v) * p * p.transpose() - pv * sm::kinematics::crossMx(p)) +
                           factor[1] * (2 * p * v.transpose() - v * p.transpose() - pv * Eigen::Matrix3d::Identity()) +
                           factor[2] * (2 * pv * p * p.transpose() - pv * p2 * Eigen::Matrix3d::Identity() - p2 * v * p.transpose());

      if(J)
    	{
    		*J = -Jpp * Jv - Jo * Ja;
    	}

      return acce;
    }

    void BSplinePose::initPoseSpline(double t0, double t1, const Eigen::Matrix4d & T_n_t0, const Eigen::Matrix4d & T_n_t1)
    {
      Eigen::VectorXd v0 = transformationToCurveValue(T_n_t0);
      Eigen::VectorXd v1 = transformationToCurveValue(T_n_t1);
      
      initSpline(t0,t1,v0,v1);
    }
    
    void BSplinePose::addPoseSegment(double tk, const Eigen::Matrix4d & T_n_tk)
    {
      Eigen::VectorXd vk = transformationToCurveValue(T_n_tk);
      
      addCurveSegment(tk, vk);
    }

    void BSplinePose::addPoseSegment2(double tk, const Eigen::Matrix4d & T_n_tk, double lambda)
    {
      Eigen::VectorXd vk = transformationToCurveValue(T_n_tk);
      
      addCurveSegment2(tk, vk, lambda);
    }


    Eigen::Matrix4d BSplinePose::curveValueToTransformation( const Eigen::VectorXd & c ) const
    {
      SM_ASSERT_EQ_DBG(Exception, c.size(), 6, "The curve value is an unexpected size!");
      Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
      T.topLeftCorner<3,3>() = rotation_->parametersToRotationMatrix(c.tail<3>());
      T.topRightCorner<3,1>() = c.head<3>();
      
      return T;
    }

    // dT/dp s.t. T(p + \delta) = Exp(dT/dp \delta)T(p)
    Eigen::Matrix4d BSplinePose::curveValueToTransformationAndJacobian( const Eigen::VectorXd & p, Eigen::MatrixXd * J ) const
    {
      SM_ASSERT_EQ_DBG(Exception, p.size(), 6, "The curve value is an unexpected size!");
      Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
      Eigen::Matrix3d S;
      T.topLeftCorner<3,3>() = rotation_->parametersToRotationMatrix(p.tail<3>(), &S);
      T.topRightCorner<3,1>() = p.head<3>();

      if(J)
      {    
        *J = Eigen::MatrixXd::Identity(6,6);
        J->topRightCorner<3,3>() = -crossMx(p.head<3>()) * S;
        J->bottomRightCorner<3,3>() = S;
      }
      
      return T;

    }

    Eigen::VectorXd BSplinePose::transformationToCurveValue( const Eigen::Matrix4d & T ) const
    {
      Eigen::VectorXd c(6);
      c.head<3>() = T.topRightCorner<3,1>();
      c.tail<3>() = rotation_->rotationMatrixToParameters(T.topLeftCorner<3,3>());
      
      return c;
    }

    void BSplinePose::initPoseSpline2(const Eigen::VectorXd & times, const Eigen::Matrix<double,6,Eigen::Dynamic> & poses, int numSegments, double lambda)
    {
      initSpline2(times, poses, numSegments, lambda);
    }

    void BSplinePose::initPoseSpline3(const Eigen::VectorXd & times, const Eigen::Matrix<double,6,Eigen::Dynamic> & poses, int numSegments, double lambda)
    {
      initSpline3(times, poses, numSegments, lambda);
    }

    void BSplinePose::initPoseSplineSparse(const Eigen::VectorXd & times, const Eigen::Matrix<double,6,Eigen::Dynamic> & poses, int numSegments, double lambda)
    {
      initSplineSparse(times, poses, numSegments, lambda);
    }
    
    void BSplinePose::initPoseSplineSparseKnots(const Eigen::VectorXd &times, const Eigen::MatrixXd &interpolationPoints, const Eigen::VectorXd knots, double lambda)
    {
    	initSplineSparseKnots(times, interpolationPoints, knots, lambda);
    }
    
    RotationalKinematics::Ptr BSplinePose::rotation() const
    {
      return rotation_;
    }
  } // namespace bsplines
