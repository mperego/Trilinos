// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_THYRAVECTOR_H
#define ROL_THYRAVECTOR_H

#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"
#include "ROL_Vector.hpp"

#include <exception>

/** \class ROL::ThyraVector
    \brief Implements the ROL::Vector interface for a Thyra Vector.
*/

namespace ROL {

template <class Real>
class ThyraVector : public Vector<Real> {
private:

  Teuchos::RCP<Thyra::VectorBase<Real> >  thyra_vec_;

  class GetEleAccelerator {
    Teuchos::RCP<const Thyra::VectorBase<Real> > vec_;

    std::vector<Teuchos::ArrayRCP<const Real> > flatVec_;

    GetEleAccelerator(); // hide default constructor

    // This function assumes Thyra vectors are just product and SPMD vectors and
    // produces a flat array of ArrayRCP objects
    std::vector<Teuchos::ArrayRCP<const Real> > buildFlatStructure(const Thyra::VectorBase<Real> & vec)
    { 
      using Teuchos::Ptr;
      using Teuchos::ptrFromRef;
      using Teuchos::ptr_dynamic_cast;

      // is an spmd vector?
      Ptr<const Thyra::SpmdVectorBase<Real> > spmd_vec = ptr_dynamic_cast<const Thyra::SpmdVectorBase<Real> >(ptrFromRef(vec));

      if(spmd_vec!=Teuchos::null) {
        std::vector<Teuchos::ArrayRCP<const Real> > flatVec(1);
        Teuchos::ArrayRCP<const Real> ptrArray;
        spmd_vec->getLocalData(ptrFromRef(ptrArray));

        flatVec[0] = ptrArray;

        return flatVec;
      }

      // it must be a product vector then
      Ptr<const Thyra::ProductVectorBase<Real> > prod_vec = ptr_dynamic_cast<const Thyra::ProductVectorBase<Real> >(ptrFromRef(vec));

      if(prod_vec!=Teuchos::null) {

        std::vector<Teuchos::ArrayRCP<const Real> > flatVec;
        for(int i=0;i<prod_vec->productSpace()->numBlocks();i++) {
          Teuchos::RCP<const Thyra::VectorBase<Real> > block = prod_vec->getVectorBlock(i);

          std::vector<Teuchos::ArrayRCP<const Real> > subVec = buildFlatStructure(*block);
          flatVec.insert(flatVec.end(),subVec.begin(),subVec.end());
        }

        return flatVec;
      }

      TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error,
          std::endl << "ROL::ThyraVector:  " <<
          "Only handling Thyra::ProductVector and Thyra::SpmdVector." <<
          std::endl);
    }

  public: 
    GetEleAccelerator(const Teuchos::RCP<const Thyra::VectorBase<Real> > & vec)
      : vec_(vec) 
    { 
      flatVec_ = buildFlatStructure(*vec);
    }

    ::Thyra::Ordinal getSize() const { return vec_->space()->dim(); }
    ::Thyra::Ordinal getLocalSize() const 
    {
      Thyra::Ordinal sum=0;
      for(std::size_t b=0;b<flatVec_.size();b++)
        sum += flatVec_[b].size(); 
      return sum;
    }

    Real getEle(::Thyra::Ordinal i) const
    { return ::Thyra::get_ele(*vec_,i); }

    Real getLocalEle(::Thyra::Ordinal i) const
    { 
      Teuchos::Ordinal sum = 0;
      std::size_t b=0;
      for(b=0;b<flatVec_.size();b++) {
        if(i<sum+flatVec_[b].size())
          break;
        sum += flatVec_[b].size();
      }

      TEUCHOS_ASSERT(b<flatVec_.size());

      return flatVec_[b][i-sum];
    }
  };

  class SetGetEleAccelerator {
    Teuchos::RCP<Thyra::VectorBase<Real> > vec_;

    SetGetEleAccelerator(); // hide default constructor

    std::vector<Teuchos::ArrayRCP<Real> > flatVec_;

    // This function assumes Thyra vectors are just product and SPMD vectors and
    // produces a flat array of ArrayRCP objects
    std::vector<Teuchos::ArrayRCP<Real> > buildFlatStructure(Thyra::VectorBase<Real> & vec)
    { 
      using Teuchos::Ptr;
      using Teuchos::ptrFromRef;
      using Teuchos::ptr_dynamic_cast;

      // is an spmd vector?
      Ptr<Thyra::SpmdVectorBase<Real> > spmd_vec = ptr_dynamic_cast<Thyra::SpmdVectorBase<Real> >(ptrFromRef(vec));

      if(spmd_vec!=Teuchos::null) {
        std::vector<Teuchos::ArrayRCP<Real> > flatVec(1);
        Teuchos::ArrayRCP<Real> ptrArray;
        spmd_vec->getNonconstLocalData(Teuchos::ptrFromRef(ptrArray));

        flatVec[0] = ptrArray;
        return flatVec;
      }

      // it must be a product vector then
      Ptr<Thyra::ProductVectorBase<Real> > prod_vec = ptr_dynamic_cast<Thyra::ProductVectorBase<Real> >(ptrFromRef(vec));

      std::vector<Teuchos::ArrayRCP<Real> > flatVec;
      for(int i=0;i<prod_vec->productSpace()->numBlocks();i++) {
        Teuchos::RCP<Thyra::VectorBase<Real> > block = prod_vec->getNonconstVectorBlock(i);

        std::vector<Teuchos::ArrayRCP<Real> > subVec = buildFlatStructure(*block);
        flatVec.insert(flatVec.end(),subVec.begin(),subVec.end());
      }

      return flatVec;
    }
  public: 
    SetGetEleAccelerator(const Teuchos::RCP<Thyra::VectorBase<Real> > & vec)
      : vec_(vec) 
    { 
      flatVec_ = buildFlatStructure(*vec);
      flatVec_ = buildFlatStructure(*vec);
    }

    ::Thyra::Ordinal getSize() const { return vec_->space()->dim(); }
    ::Thyra::Ordinal getLocalSize() const 
    {
      Thyra::Ordinal sum=0;
      for(std::size_t b=0;b<flatVec_.size();b++)
        sum += flatVec_[b].size(); 
      return sum;
    }

    void setEle(::Thyra::Ordinal i,Real v) 
    { ::Thyra::set_ele(i,v,vec_.ptr()); }

    void setLocalEle(::Thyra::Ordinal i,Real v) const
    { 
      Teuchos::Ordinal sum = 0;
      std::size_t b=0;
      for(b=0;b<flatVec_.size();b++) {
        if(i<sum+flatVec_[b].size())
          break;
        sum += flatVec_[b].size();
      }

      TEUCHOS_ASSERT(b<flatVec_.size());

      flatVec_[b][i-sum] = v;
    }

    Real getEle(::Thyra::Ordinal i) const
    { return ::Thyra::get_ele(*vec_,i); }

    Real getLocalEle(::Thyra::Ordinal i) const
    { 
      Teuchos::Ordinal sum = 0;
      std::size_t b=0;
      for(b=0;b<flatVec_.size();b++) {
        if(i<sum+flatVec_[b].size())
          break;
        sum += flatVec_[b].size();
      }

     
      std::stringstream ss;
      ss << "Block identifier b= " << b << " is too large for i=" << i << " on array with " << getLocalSize() <<
            " and " << flatVec_.size() << " blocks.";
      ROL_TEST_FOR_EXCEPTION(b>=flatVec_.size(),std::logic_error, ss.str());

      return flatVec_[b][i-sum];
    }
  };

public:
  ~ThyraVector() {}

  ThyraVector(const Teuchos::RCP<Thyra::VectorBase<Real> > & thyra_vec) : thyra_vec_(thyra_vec) {}

  /** \brief Compute \f$y \leftarrow x + y\f$ where \f$y = \mbox{*this}\f$.
  */
  void plus( const Vector<Real> &x ) {
    const ThyraVector &ex = dynamic_cast<const ThyraVector&>(x);
    ::Thyra::Vp_V( thyra_vec_.ptr(), *ex.getVector());
  }

  /** \brief Compute \f$y \leftarrow \alpha y\f$ where \f$y = \mbox{*this}\f$.
  */
  void scale( const Real alpha ) { 
    ::Thyra::scale(alpha, thyra_vec_.ptr());
  }

  /** \brief Returns \f$ \langle y,x \rangle \f$ where \f$y = \mbox{*this}\f$.
  */
  Real dot( const Vector<Real> &x ) const {
    const ThyraVector &ex = dynamic_cast<const ThyraVector&>(x);
    return ::Thyra::dot<Real>(*thyra_vec_, *ex.thyra_vec_);
  }

  /** \brief Apply \f$\mathtt{*this}\f$ to a dual vector.  This is equivalent
             to the call \f$\mathtt{this->dot(x.dual())}\f$.
  */
  Real apply( const Vector<Real> &x ) const {
    const ThyraVector &ex = dynamic_cast<const ThyraVector&>(x);
    return ::Thyra::dot<Real>(*thyra_vec_, *ex.thyra_vec_);
  }

  /** \brief Returns \f$ \| y \| \f$ where \f$y = \mbox{*this}\f$.
  */
  Real norm() const {
    // return ::Thyra::norm_2<Real>(*thyra_vec_);
    return std::sqrt(dot(*this));
  } 

  /** \brief Clone to make a new (uninitialized) vector.
  */
  Teuchos::RCP<Vector<Real> > clone() const{
    Teuchos::RCP<Thyra::VectorBase<Real> > tv = thyra_vec_->clone_v();
    return Teuchos::rcp( new ThyraVector(tv) ); 
  }

  /** \brief Compute \f$y \leftarrow \alpha x + y\f$ where \f$y = \mbox{*this}\f$.
  */
  void axpy( const Real alpha, const Vector<Real> &x ) {
    const ThyraVector &ex = dynamic_cast<const ThyraVector&>(x);
    ::Thyra::Vp_StV( thyra_vec_.ptr(), alpha, *ex.getVector());
  }

  /**  \brief Set to zero vector.
  */
  void zero() {
    ::Thyra::put_scalar(0.0, thyra_vec_.ptr());
  }

  /**  \brief Set all entries of the vector to alpha.
    */
  void setScalar(const Real C) {
      ::Thyra::put_scalar(C, thyra_vec_.ptr());
    }

  /**  \brief Set entries of the vector to uniform random between l and u.
    */
  void randomize(const Real l=0.0, const Real u=1.0) {
      ::Thyra::randomize(l, u, thyra_vec_.ptr());
    }

  /**  \brief Set all entries of the vector to alpha.
    */
  void putScalar(Real alpha) {
      ::Thyra::put_scalar(alpha, thyra_vec_.ptr());
    }

  /**  \brief const get of the Thyra vector.
    */
  Teuchos::RCP<const Thyra::VectorBase<Real> > getVector() const {
    return thyra_vec_;
  }

  Teuchos::RCP<const Thyra::MultiVectorBase<Real> > getMultiVector() const {
    return thyra_vec_;
  }

  /**  \brief nonconst get of the Thyra vector.
    */
  Teuchos::RCP<Thyra::VectorBase<Real> > getVector()  {
    return thyra_vec_;
  }

  Teuchos::RCP<Thyra::MultiVectorBase<Real> > getMultiVector() {
    return thyra_vec_;
  }

  /** \brief Return i-th basis vector.
    */
  Teuchos::RCP<Vector<Real> > basis( const int i ) const {
    Teuchos::RCP<Vector<Real> > e = clone();
    Teuchos::RCP<Thyra::VectorBase<Real> > basisThyraVec = (Teuchos::rcp_static_cast<ThyraVector>(e))->getVector();
    ::Thyra::put_scalar(0.0, basisThyraVec.ptr());
    ::Thyra::set_ele(i,1.0, basisThyraVec.ptr());
    return e;
  }

  /** \brief Return dimension of the vector space.
    */
  int dimension() const {
    return thyra_vec_->space()->dim();
  }

  /** \brief Set \f$y \leftarrow x\f$ where \f$y = \mathtt{*this}\f$.
    */
  void set(const Vector<Real> &x ) {
    const ThyraVector &ex = dynamic_cast<const ThyraVector&>(x);
    ::Thyra::copy( *ex.getVector(), thyra_vec_.ptr() );
  }

  void applyUnary( const Elementwise::UnaryFunction<Real> &f ) {
//*
    SetGetEleAccelerator thisAccel(thyra_vec_);

    for(::Thyra::Ordinal i=0;i<thisAccel.getLocalSize();i++) {
      Real val  = thisAccel.getLocalEle(i); 

      thisAccel.setLocalEle(i,f.apply(val));
    }
/*/
    //int map[306] = {0, 1, 52, 51, 2, 53, 3, 54, 4, 55, 5, 56, 6, 57, 7, 58, 8, 59, 9, 60, 10, 61, 11, 62, 12, 63, 13, 64, 14, 65, 15, 66, 16, 67, 17, 68, 18, 69, 19, 70, 20, 71, 21, 72, 22, 73, 23, 74, 24, 75, 25, 76, 26, 77, 27, 78, 28, 79, 29, 80, 30, 81, 31, 82, 32, 83, 33, 84, 34, 85, 35, 86, 36, 87, 37, 88, 38, 89, 39, 90, 40, 91, 41, 92, 42, 93, 43, 94, 44, 95, 45, 96, 46, 97, 47, 98, 48, 99, 49, 100, 50, 101, 103, 102, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 154, 153, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 205, 204, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 256, 255, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305};
    int map[306] = {0, 1, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 3, 2, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97, 99, 101, 103, 102, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 154, 153, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 205, 204, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 256, 255, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305};
    if(thyra_vec_->space()->dim()==51){
      for(::Thyra::Ordinal i=0;i<51;i++)
        map[i]=i;
    }

    for(::Thyra::Ordinal i=0;i<51;i++) {
      Real val = ::Thyra::get_ele(*thyra_vec_,map[i]);
      ::Thyra::set_ele(map[i],f.apply(val),thyra_vec_.ptr());
    }
    for(::Thyra::Ordinal i=51;i<thyra_vec_->space()->dim();i++) {
      Real val = ::Thyra::get_ele(*thyra_vec_,map[i%51]);
      ::Thyra::set_ele(map[i],val,thyra_vec_.ptr());
    }
//*/
  }

  void applyBinary( const Elementwise::BinaryFunction<Real> &f, const Vector<Real> &x ) {
    const ThyraVector &ex = dynamic_cast<const ThyraVector&>(x);
    Teuchos::RCP< const Thyra::VectorBase<Real> > xp = ex.getVector();

    SetGetEleAccelerator thisAccel(thyra_vec_);
    GetEleAccelerator xpAccel(xp);

    TEUCHOS_ASSERT(thisAccel.getLocalSize()==xpAccel.getLocalSize());
    for(::Thyra::Ordinal i=0;i<thisAccel.getLocalSize();i++) {
      Real val  = thisAccel.getLocalEle(i); 
      Real xval = xpAccel.getLocalEle(i); 

      thisAccel.setLocalEle(i,f.apply(val,xval));
    }
/*
    TEUCHOS_ASSERT(thisAccel.getSize()==thyra_vec_->space()->dim());
    TEUCHOS_ASSERT(xpAccel.getSize()==xp->space()->dim());

    for(::Thyra::Ordinal i=0;i<thyra_vec_->space()->dim();i++) {
      Real val  = thisAccel.getEle(i); // ::Thyra::get_ele(*thyra_vec_,i);
      Real xval = xpAccel.getEle(i); // ::Thyra::get_ele(*xp,i);

      // ::Thyra::set_ele(i,f.apply(val,xval),thyra_vec_.ptr());
      thisAccel.setEle(i,f.apply(val,xval));
    }
*/
  }

  Real reduce( const Elementwise::ReductionOp<Real> &r ) const {
    Real result = r.initialValue();
    
    for(::Thyra::Ordinal i=0;i<thyra_vec_->space()->dim();i++) {
      r.reduce(::Thyra::get_ele(*thyra_vec_,i),result);  
    }
    return result;
  } 

}; // class ThyraVector

} // namespace ROL

#endif

