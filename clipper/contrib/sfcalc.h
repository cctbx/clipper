/*!  \file sfcalc.h
  Header file for sample structure factor calculation impelementation
  \ingroup g_sfcalc
*/
//C Copyright (C) 2000-2004 Kevin Cowtan and University of York
//L
//L  This library is free software and is distributed under the terms and
//L  conditions of the CCP4 licence agreement as `Part 0' (Annex 2)
//L  software, which is version 2.1 of the GNU Lesser General Public
//L  Licence (LGPL) with the following additional clause:
//L
//L     `You may also combine or link a "work that uses the Library" to
//L     produce a work containing portions of the Library, and distribute
//L     that work under terms of your choice, provided that you give
//L     prominent notice with each copy of the work that the specified
//L     version of the Library is used in it, and that you include or
//L     provide public access to the complete corresponding
//L     machine-readable source code for the Library including whatever
//L     changes were used in the work. (i.e. If you make changes to the
//L     Library you must distribute those, but you do not need to
//L     distribute source or object code to those portions of the work
//L     not covered by this licence.)'
//L
//L  Note that this clause grants an additional right and does not impose
//L  any additional restriction, and so does not affect compatibility
//L  with the GNU General Public Licence (GPL). If you wish to negotiate
//L  other terms, please contact the maintainer.
//L
//L  You can redistribute it and/or modify the library under the terms of
//L  the GNU Lesser General Public License as published by the Free Software
//L  Foundation; either version 2.1 of the License, or (at your option) any
//L  later version.
//L
//L  This library is distributed in the hope that it will be useful, but
//L  WITHOUT ANY WARRANTY; without even the implied warranty of
//L  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//L  Lesser General Public License for more details.
//L
//L  You should have received a copy of the CCP4 licence and/or GNU
//L  Lesser General Public License along with this library; if not, write
//L  to the CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//L  The GNU Lesser General Public can also be obtained by writing to the
//L  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//L  MA 02111-1307 USA


#ifndef CLIPPER_SFCALC
#define CLIPPER_SFCALC


#include "function_object_bases.h"


namespace clipper {


  //! Isotropic structure factor calculation by slow summation
  /*! \ingroup g_sfcalc */
  template<class T> class SFcalc_iso_sum : public SFcalc_base<T> {
  public:
    //! constructor
    SFcalc_iso_sum() {}
    //! constructor: shorthand for constructor+operator
    SFcalc_iso_sum( HKL_data<datatypes::F_phi<T> >& fphidata, const Atom_list& atoms ) { (*this)( fphidata, atoms ); }
    bool operator() ( HKL_data<datatypes::F_phi<T> >& fphidata, const Atom_list& atoms ) const;
  };


  //! Anisotropic structure factor calculation by slow summation
  /*! \ingroup g_sfcalc */
  template<class T> class SFcalc_aniso_sum : public SFcalc_base<T> {
  public:
    //! constructor
    SFcalc_aniso_sum() {}
    //! constructor: shorthand for constructor+operator
    SFcalc_aniso_sum( HKL_data<datatypes::F_phi<T> >& fphidata, const Atom_list& atoms ) { (*this)( fphidata, atoms ); }
    bool operator() ( HKL_data<datatypes::F_phi<T> >& fphidata, const Atom_list& atoms ) const;
  };


  //! Isotropic structure factor calculation by fast Fourier
  /*! \ingroup g_sfcalc */
  template<class T> class SFcalc_iso_fft : public SFcalc_base<T> {
  public:
    //! constructor
    /*! \param rate Shannon rate (oversampling) of the FFT grid.
        \param uadd Additional U for smoothing atoms. */
    SFcalc_iso_fft( const ftype radius = 2.5, const ftype rate = 1.5, const ftype uadd = 0.0 ) : radius_(radius), rate_(rate), uadd_(uadd) {}
    //! constructor: shorthand for constructor+operator
    SFcalc_iso_fft( HKL_data<datatypes::F_phi<T> >& fphidata, const Atom_list& atoms, const ftype radius = 2.5, const ftype rate = 1.5, const ftype uadd = 0.0 ) : radius_(radius), rate_(rate), uadd_(uadd) { (*this)( fphidata, atoms ); }
    bool operator() ( HKL_data<datatypes::F_phi<T> >& fphidata, const Atom_list& atoms ) const;
  private:
    const ftype radius_, rate_, uadd_;
  };


  //! Anisotropic structure factor calculation by fast Fourier
  /*! \ingroup g_sfcalc */
  template<class T> class SFcalc_aniso_fft : public SFcalc_base<T> {
  public:
    //! constructor
    /*! \param rate Shannon rate (oversampling) of the FFT grid.
        \param uadd Additional U for smoothing atoms. */
    SFcalc_aniso_fft( const ftype radius = 2.5, const ftype rate = 1.5, const ftype uadd = 0.0 ) : radius_(radius), rate_(rate), uadd_(uadd) {}
    //! constructor: shorthand for constructor+operator
    SFcalc_aniso_fft( HKL_data<datatypes::F_phi<T> >& fphidata, const Atom_list& atoms, const ftype radius = 2.5, const ftype rate = 1.5, const ftype uadd = 0.0 ) : radius_(radius), rate_(rate), uadd_(uadd) { (*this)( fphidata, atoms ); }
    bool operator() ( HKL_data<datatypes::F_phi<T> >& fphidata, const Atom_list& atoms ) const;
  private:
    const ftype radius_, rate_, uadd_;
  };


} // namespace clipper

#endif
