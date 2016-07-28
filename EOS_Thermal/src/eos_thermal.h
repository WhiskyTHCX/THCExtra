/*! \file eos_thermal.h
\brief Defines generic API for thermal EOS
\author Wolfgang Kastaun

This contains classes representing a generic thermal EOS, where the independent
quantities are rest mass density, specific energy, and electron fraction.

The class eos_thermal provides all functions defining the EOS, but does
not implement any particular EOS by itself. Instead it holds a pointer to
to an implementation to which it forwards the calls.
The user does not need to know, and in fact cannot know, the actual type of the EOS.
Implementations are derived from the abstract polymorphic base class eos_thermal_impl.
The latter defines the interface one has to implement when writing a new EOS.
*/

#ifndef EOSTHERMAL_H
#define EOSTHERMAL_H
#include <string>
#include <boost/shared_ptr.hpp>

namespace whizza {

typedef double pz_real;

/// Class representing error conditions in EOS calls.
struct eos_thermal_status {
  bool failed;          ///< Set to true if parameters are out of range/NAN/INF
  std::string err_msg;  ///< Error description in case of failure, else undefined.
  /// Default constructor: Set to no failure.
  eos_thermal_status() : failed(false) {}
  /// Set fail flag and error message.
  void fail(std::string msg) {failed=true; err_msg=msg;}
};


///Class representing a range
struct eos_thermal_range {
  pz_real min;  ///< Minimum
  pz_real max;  ///< Maximum
  ///Default constructor: empty range.
  eos_thermal_range() : min(0), max(0) {}
  ///Construct from minimum and maximum
  eos_thermal_range(pz_real min_, pz_real max_) : min(min_), max(max_) {}
  ///Check if value is contained in [min,max]. False for NAN or INF.
  bool contains(pz_real x) const {return (x>=min) && (x<=max);}
};


/// Virtual abstract base class defining EOS implementation API.
/**
All implementations of thermal EOS are derived from this class.
This class also defines the validity ranges of the EOS, which
have to be provided by the implementations.
The EOS functions of the implementations like press_from_valid_rho_eps_ye
can assume that the ranges of the parameters are valid, since they are
protected and never called directly. Instead, this class provides
non-virtual functions, e.g. press_from_rho_eps_ye, which check the ranges
and then call the virtual functions of the implementation.
The public functions of this class all take a reference to a structure
eos_thermal_impl::status for error handling. In case of range errors,
the failed flag is set, the results are explicitly set to NAN, and the
status::err_msg field contains a string with a description of the error.
Otherwise, the failed flag is cleared (err_msg is left untouched for
performance reasons, string operations are slow).
Implementations should always succeed when fed with parameters in the
valid range. A failure to do so is an abnormal condition and should be
handled by throwing an exception or program termination.
Exceptions from the implementation will not be caught by this class.
**/
class eos_thermal_impl {
  public:
  typedef eos_thermal_status status;
  typedef eos_thermal_range range;

  private:

  range rgrho;    ///< Valid range for density \f$ \rho \f$
  range rgye;     ///< Valid range for electron fraction \f$ Y_e \f$
  range rgtemp;   ///< Valid range for temperature \f$ T \f$

  public:

  eos_thermal_impl() {}
  virtual ~eos_thermal_impl();

  /// Compute pressure
  /**
  \returns Pressure in the range \f$ 0 \le P < \rho (1+\epsilon) \f$
  **/
  pz_real press_from_rho_eps_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real eps,     ///<Specific internal energy \f$ \epsilon \f$
    const pz_real ye,      ///<Electron fraction \f$ Y_e \f$
    status&       stat     ///<Error report
  ) const;

  /// Compute speed of sound.
  /**
  \returns Soundspeed in the range \f$  0 \le c_s < 1  \f$
  **/
  pz_real csnd_from_rho_eps_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real eps,     ///<Specific internal energy \f$ \epsilon \f$
    const pz_real ye,      ///<Electron fraction \f$ Y_e \f$
    status&       stat     ///<Error report
  ) const;

  /// Compute temperature
  /**
  \returns Temperature in MeV.
  **/
  pz_real temp_from_rho_eps_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real eps,     ///<Specific internal energy \f$ \epsilon \f$
    const pz_real ye,      ///<Electron fraction \f$ Y_e \f$
    status&       stat     ///<Error report
  ) const;

  /// Compute pressure and derivatives.
  void press_derivs_from_rho_eps_ye(
    pz_real& press,        ///<Pressure \f$ P \f$
    pz_real& dpdrho,       ///<Partial derivative \f$ \frac{\partial P}{\partial \rho} \f$
    pz_real& dpdeps,       ///<Partial derivative \f$ \frac{\partial P}{\partial \epsilon} \f$
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real eps,     ///<Specific internal energy \f$ \epsilon \f$
    const pz_real ye,      ///<Electron fraction \f$ Y_e \f$
    status&       stat     ///<Error report
  ) const;

  /// Compute pressure and soundspeed.
  void press_csnd_from_rho_eps_ye(
    pz_real& press,        ///<Pressure \f$ P \f$
    pz_real& csnd,         ///<Soundspeed \f$ c_s \f$
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real eps,     ///<Specific internal energy \f$ \epsilon \f$
    const pz_real ye,      ///<Electron fraction \f$ Y_e \f$
    status&       stat     ///<Error report
  ) const;

  /// Compute specific entropy
  pz_real entropy_from_rho_eps_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real eps,     ///<Specific internal energy \f$ \epsilon \f$
    const pz_real ye,      ///<Electron fraction \f$ Y_e \f$
    status&       stat     ///<Error report
  ) const;


  /// Compute specific energy
  pz_real eps_from_rho_temp_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real temp,    ///<Temperature \f$ T \f$ in MeV
    const pz_real ye,      ///<Electron fraction \f$ Y_e \f$
    status&       stat     ///<Error report
  ) const;

  /// Compute pressure
  pz_real press_from_rho_temp_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real temp,    ///<Temperature
    const pz_real ye,      ///<Electron fraction \f$ Y_e \f$
    status&       stat     ///<Error report
  ) const;

  /// Compute speed of sound.
  pz_real csnd_from_rho_temp_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real temp,    ///<Temperature
    const pz_real ye,      ///<Electron fraction \f$ Y_e \f$
    status&       stat     ///<Error report
  ) const;

  /// Compute pressure and soundspeed.
  void press_csnd_from_rho_temp_ye(
    pz_real& press,        ///<Pressure \f$ P \f$
    pz_real& csnd,         ///<Soundspeed \f$ c_s \f$
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real temp,    ///<Temperature
    const pz_real ye,      ///<Electron fraction \f$ Y_e \f$
    status&       stat     ///<Error report
  ) const;

  /// Compute specific entropy
  pz_real entropy_from_rho_temp_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real temp,    ///<Temperature \f$ T \f$ in MeV
    const pz_real ye,      ///<Electron fraction \f$ Y_e \f$
    status&       stat     ///<Error report
  ) const;

  /// Valid range for specifc energy at given density and electron fraction.
  range range_eps(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const;

  /// Valid range for density
  const range& range_rho() const {return rgrho;}

  /// Valid range for electron fraction
  const range& range_ye() const {return rgye;}

  /// Valid range for temperature
  const range& range_temp() const {return rgtemp;}

  /// Check if density is in the validity range of the EOS.
  bool is_rho_valid(pz_real rho) const {return range_rho().contains(rho);}

  /// Check if electron fraction is in the validity range of the EOS.
  bool is_ye_valid(pz_real ye) const {return range_ye().contains(ye);}

  /// Check if temperature is in the validity range of the EOS.
  bool is_temp_valid(pz_real temp) const {return range_temp().contains(temp);}

  /// Check if specific energy is in the valid range at given density and electron fraction.
  bool is_eps_valid(pz_real rho, pz_real eps, pz_real ye) const {return range_eps(rho, ye).contains(eps);}

  /// Check if density, specific energy, and electron fraction are in the valid range.
  bool is_rho_eps_ye_valid(pz_real rho, pz_real eps, pz_real ye) const
  {
    return (is_rho_valid(rho) && is_ye_valid(ye) && is_eps_valid(rho, eps, ye));
  }

  /// Check if density, temperature, and electron fraction are in the valid range.
  bool is_rho_temp_ye_valid(pz_real rho, pz_real temp, pz_real ye)  const
  {
    return (is_rho_valid(rho) && is_ye_valid(ye) && is_temp_valid(temp));
  }

  /// Check for range error, set status report, provide detailed error message in case of error.
  void check_rho_eps_ye(const pz_real rho, const pz_real eps, const pz_real ye, status& stat) const;

  /// Check for range error, set status report, provide detailed error message in case of error.
	void check_rho_temp_ye(const pz_real rho, const pz_real temp, const pz_real ye, status& stat) const;

  /// Helper function: returns a NAN
  static pz_real nan();


  /// Compute pressure for valid input. To be implemented by derived class.
  /** Implementations should always guarantee
      \f$ 0 \le P < \rho (1+\epsilon) \f$
  **/
  virtual pz_real press_from_valid_rho_eps_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real eps,     ///<Specific internal energy \f$ \epsilon \f$
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const =0;

  /// Compute sound speed for valid input. To be implemented by derived class.
  /** Implementations should always guarantee
      \f$  0 \le c_s < 1  \f$
  **/
  virtual pz_real csnd_from_valid_rho_eps_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real eps,     ///<Specific internal energy \f$ \epsilon \f$
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const =0;

  /// Compute temperature for valid input. To be implemented by derived class.
  /**
  \returns Temperature in MeV
  **/
  virtual pz_real temp_from_valid_rho_eps_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real eps,     ///<Specific internal energy \f$ \epsilon \f$
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const =0;

  /// Compute pressure and derivatives for valid input. To be implemented by derived class.
  virtual void press_derivs_from_valid_rho_eps_ye(
    pz_real& press,        ///<Pressure \f$ P \f$
    pz_real& dpdrho,       ///<Partial derivative \f$ \frac{\partial P}{\partial \rho} \f$
    pz_real& dpdeps,       ///<Partial derivative \f$ \frac{\partial P}{\partial \epsilon} \f$
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real eps,     ///<Specific internal energy \f$ \epsilon \f$
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const =0;

  /// Compute specific entropy for valid input. To be implemented by derived class.
  virtual pz_real entropy_from_valid_rho_eps_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real eps,     ///<Specific internal energy \f$ \epsilon \f$
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const =0;


  /// Compute pressure for valid input.
  virtual pz_real press_from_valid_rho_temp_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real temp,     ///<Temperature
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const;

  /// Compute sound speed for valid input.
  virtual pz_real csnd_from_valid_rho_temp_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real temp,     ///<Temperature
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const;


  /// Compute specific entropy for valid input.
  virtual pz_real entropy_from_valid_rho_temp_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real temp,    ///<Temperature \f$ T \f$
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const;

  /// Compute specific energy for valid input. To be implemented by derived class.
  virtual pz_real eps_from_valid_rho_temp_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real temp,    ///<Temperature \f$ T \f$
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const =0;

  /// Compute the valid range of specific energy for given density and electron fraction.
  virtual range range_eps_from_valid_rho_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const =0;

  protected:

  /// Set the density range. Has to be called in the constructor of any implementation.
  void set_range_rho(const range& r) {rgrho = r;}
  /// Set the electron fraction range. Has to be called in the constructor of any implementation.
  void set_range_ye(const range& r) {rgye = r;}
  /// Set the temperature range. Has to be called in the constructor of any implementation.
  void set_range_temp(const range& r) {rgtemp = r;}

};


/// Class representing a generic thermal EOS
/**
This is a wrapper holding a pointer to the actual implementation of the EOS.
This pointer is a reference counted shared pointer.
Objects of this class can be copied or assigned without worries, since
only the pointer to the implementation is copied, and the implementation is deleted
once the last object pointing to it is destroyed.
Further, this wrapper only calls const methods of the implementation it holds,
which remains unchanged until deletion.

This wrapper provides the same public interface as eos_thermal_impl, see there
for documentation.
**/
class eos_thermal {
  boost::shared_ptr<eos_thermal_impl> pimpl;  ///< Shared pointer to implementation
  public:
  typedef eos_thermal_impl::range range;
  typedef eos_thermal_impl::status status;
  const eos_thermal_impl& implementation() const;

///Default constructor
/**
Creates uninitialized EOS. Any attemt to use resulting object throws an exception.
One can use it after copying another EOS to it.
*/
  eos_thermal() {}

///Constructor from pointer to implementation.
/**
The object takes ownership of the pointer, do not delete the pointer.
*/
  eos_thermal(
    eos_thermal_impl* eos   ///< Pointer to implementation
  ) : pimpl(eos) {}

///Copy constructor.
/**
This does not copy the implementation itself, e.g. a large tabulated EOS
has to be kept in memory only once.
**/
  eos_thermal(const eos_thermal& that) : pimpl(that.pimpl) {}

///Assignment operator.
  eos_thermal& operator=(const eos_thermal& that)
  {
    pimpl=that.pimpl; return *this;
  }

  pz_real press_from_rho_eps_ye(const pz_real rho,const pz_real eps,
     const pz_real ye, status& stat) const
  {
    return implementation().press_from_rho_eps_ye(rho, eps, ye, stat);
  }

  pz_real csnd_from_rho_eps_ye(const pz_real rho, const pz_real eps,
      const pz_real ye, status& stat) const
  {
    return implementation().csnd_from_rho_eps_ye(rho, eps, ye, stat);
  }

  pz_real temp_from_rho_eps_ye(const pz_real rho, const pz_real eps,
      const pz_real ye, status& stat) const
  {
    return implementation().temp_from_rho_eps_ye(rho, eps, ye, stat);
  }

  void press_derivs_from_rho_eps_ye(pz_real& press, pz_real& dpdrho,
            pz_real& dpdeps, const pz_real rho, const pz_real eps,
            const pz_real ye, status& stat) const
  {
    return implementation().press_derivs_from_rho_eps_ye(press,
                                 dpdrho, dpdeps, rho, eps, ye, stat);
  }

  pz_real entropy_from_rho_temp_ye(const pz_real rho, const pz_real temp,
      const pz_real ye, status& stat) const
  {
    return implementation().entropy_from_rho_temp_ye(rho, temp, ye, stat);
  }

  pz_real entropy_from_rho_eps_ye(const pz_real rho, const pz_real eps,
      const pz_real ye, status& stat) const
  {
    return implementation().entropy_from_rho_eps_ye(rho, eps, ye, stat);
  }

  void press_csnd_from_rho_eps_ye(pz_real& press, pz_real& csnd,
    const pz_real rho, const pz_real eps, const pz_real ye,
    status& stat) const
  {
    return implementation().press_csnd_from_rho_eps_ye(press, csnd, rho,
                                 eps, ye, stat);
  }

  pz_real eps_from_rho_temp_ye(const pz_real rho, const pz_real temp,
      const pz_real ye, status& stat) const
  {
    return implementation().eps_from_rho_temp_ye(rho, temp, ye, stat);
  }

  range range_eps(const pz_real rho, const pz_real ye) const
  {
    return implementation().range_eps(rho, ye);
  }

  const range& range_rho() const
  {
    return implementation().range_rho();
  }

  const range& range_ye() const
  {
    return implementation().range_ye();
  }

  const range& range_temp() const
  {
    return implementation().range_temp();
  }

  bool is_rho_valid(pz_real rho) const
  {
    return implementation().is_rho_valid(rho);
  }

  bool is_ye_valid(pz_real ye) const
  {
    return implementation().is_ye_valid(ye);
  }

  bool is_temp_valid(pz_real temp) const
  {
    return implementation().is_temp_valid(temp);
  }

  bool is_eps_valid(pz_real rho, pz_real eps, pz_real ye) const
  {
    return implementation().is_eps_valid(rho, eps, ye);
  }

  bool is_rho_eps_ye_valid(pz_real rho, pz_real eps, pz_real ye) const
  {
    return implementation().is_rho_eps_ye_valid(rho, eps, ye);
  }

  bool is_rho_temp_ye_valid(pz_real rho, pz_real temp, pz_real ye)  const
  {
    return implementation().is_rho_temp_ye_valid(rho, temp, ye);
  }

  void check_rho_eps_ye(const pz_real rho, const pz_real eps,
                        const pz_real ye, status& stat) const
  {
    implementation().check_rho_eps_ye(rho, eps, ye, stat);
  }

	void check_rho_temp_ye(const pz_real rho, const pz_real temp,
                         const pz_real ye, status& stat) const
  {
    implementation().check_rho_temp_ye(rho, temp, ye, stat);
  }


  static pz_real nan() {return eos_thermal_impl::nan();}
};

}

#endif

