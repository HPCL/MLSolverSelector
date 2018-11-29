#ifndef DEFAULTPETSCFEATURE_SET
#define DEFAULTPETSCFEATURE_SET

#ifdef WITH_PETSCUI 

#include "typedefs.h"
#include "petscksp.h"
#include "petsc/private/kspimpl.h"
#include "petscdmda.h"
#include <thread>
#include <future>

namespace SolverSelecter
{


class RaplMeasurements {
  public:
  // Reset the class so it is ready to start measuring
  std::chrono::nanoseconds m_interval;
  std::promise<void> *exitSignal = nullptr;
  std::thread *measurementThread = nullptr;

  //std::unique_ptr<Rapl> rapl_ptr; 

  RaplMeasurements(std::chrono::nanoseconds interval) : m_interval(interval) {
  
  }

  virtual void reset() {
  
  }

  virtual void startPowerMonitoring() {
       exitSignal = new std::promise<void>();
       std::future<void> futureObj = exitSignal->get_future();
       measurementThread = new std::thread(&RaplMeasurements::PowerMeasurementThreadFunction, this, std::move(futureObj));       
  
       //rapl_ptr.reset( new Rapl() );   
  }

  virtual void stopPowerMonitoring() {
       exitSignal->set_value();
       measurementThread->join();
       delete exitSignal;   
       delete measurementThread;    
  };
  
  virtual void PowerMeasurementThreadFunction(std::future<void> futureObj) {
      std::cout << "Thread Start" << std::endl;
	    while (futureObj.wait_for(m_interval) == std::future_status::timeout)
	    {
		      std::cout << " \t THREAD : Doing Some Work" << std::endl;
          //rapl->sample();
      }
	    std::cout << "Thread End" << std::endl;
  
  }

};


class PetscTestingSpace
{
public:
  

    PetscTestingSpace( );

    virtual void
    extract_features(KSP &ksp,
                     std::map<std::string, double> &fmap,
                     int edge,
                     int internal,
                     bool matvec ) = 0;

    virtual void
    start_measurements(KSP &ksp,
                       Vec &x,
                       Vec &b ) = 0;

    virtual void
    stop_measurements(KSP &ksp,
                      Vec &x,
                      Vec &b,
                      std::map<std::string, double > &mmap) = 0;

    virtual void
    classify_measurements(std::map<std::string, double> &cmap ) = 0;

};


class DefaultPetscTestingSpace : public PetscTestingSpace
{
public:

    std::chrono::time_point<std::chrono::high_resolution_clock> time_start;
    PetscLogDouble space;
    RaplMeasurements raplMeasurement;    

    bool measuring = false;

    DefaultPetscTestingSpace(std::chrono::nanoseconds interval);

    void
    extract_features( KSP &ksp,
                      std::map<std::string, double> &fmap,
                      int edge,
                      int interior,
                      bool matvecs ) override;
    void
    start_measurements(KSP &ksp,
                       Vec &x,
                       Vec &b ) override;

    void
    stop_measurements(KSP &ksp,
                      Vec &x,
                      Vec &b,
                      std::map<std::string, double> &mmap ) override;

    void
    classify_measurements(std::map<std::string, double> &cmap ) override;

    void powerMeasurementTask(); 


};

}
#endif
#endif
