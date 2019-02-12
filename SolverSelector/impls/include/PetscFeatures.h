#ifndef DEFAULTPETSCFEATURE_SET
#define DEFAULTPETSCFEATURE_SET

#ifdef WITH_PETSCUI 

#include "typedefs.h"
#include "petscksp.h"
#include "petsc/private/kspimpl.h"
#include "petscdmda.h"
#include "Rapl.h"
#include "papi.h"
#include "papi_test.h"
#include <thread>
#include <future>


#define MAX_RAPL_EVENTS 64

namespace SolverSelecter
{


class MemClass {
  public:

  int num_events=6;
  long_long *values;
  int Events[6] = { PAPI_TOT_INS, PAPI_TOT_CYC, PAPI_L1_TCM, PAPI_L1_DCM, PAPI_L2_DCA}, EventSet = PAPI_NULL;
  int retval;


   virtual void startMemMonitoring() 
   {

	/* PAPI Initialization */
     retval = PAPI_library_init( PAPI_VER_CURRENT );
     if ( retval != PAPI_VER_CURRENT ) {
  	printf("PAPI_library_init failed %d \n",retval);
     }

    values=(long_long *)calloc(num_events,sizeof(long_long));
    retval = PAPI_create_eventset(&EventSet);


    retval = PAPI_add_events(EventSet, Events, num_events);
    retval = PAPI_start(EventSet);

   }

    virtual void stopMemMonitoring( std::map<std::string, double> &mmap) 
    {
    retval = PAPI_stop(EventSet, values);
    printf("Events %d %d %d %d %d %d \n", Events[0], Events[1], Events[2], Events[3], Events[4], Events[5]);
    printf("PAPI_TOT_INS %.4llu PAPI_TOT_CYC %.4llu PAPI_L1_TCM %.4llu PAPI_L1_DCM %.4llu PAPI_L2_DCM %.4llu PAPI_L1_DCA %.4llu\n",values[0], values[1], values[2], values[3], values[4], values[5]);
    printf("------------------------------\n");
    printf("Hardware Counters   : %d\n", PAPI_num_counters());
    printf("CPI                 : %.4lf\n", values[1] * 1.0 / values[0]);
   // printf("L1 miss rate        : %.4lf\n", values[2] * 1.0 / values[3]);
   // printf("L2 miss rate        : %.4lf\n", values[4] * 1.0 / values[2]);
    }

};

class CacheClass {
  public:

  int num_event=6;
  long_long *values;
  int Events[6] = {  PAPI_L1_DCM, PAPI_L2_DCM}, EventSet = PAPI_NULL;
  int retval;

   virtual void startMemMonitoring() 
   {
     	/* PAPI Initialization */
     retval = PAPI_library_init( PAPI_VER_CURRENT );
     if ( retval != PAPI_VER_CURRENT ) {
  	printf("PAPI_library_init failed %d \n",retval);
     }


    values=(long_long *)calloc(num_event,sizeof(long_long));
    retval = PAPI_create_eventset(&EventSet);

    retval = PAPI_add_events(EventSet, Events, num_event);
    retval = PAPI_start(EventSet);

   }

    virtual void stopMemMonitoring( std::map<std::string, double> &mmap) 
    {
    retval = PAPI_stop(EventSet, values);
    printf("PAPI_L1_DCM %.4llu PAPI_L2_DCM %.4llu  \n",values[0], values[1]);
    printf("------------------------------\n");
    printf("L2 cache data hit rate   :  %.4lf\n", 1.0 - (values[1]/ values[0]));
    //printf("L1 cache data hit rate    : %.4lf\n", 1.0 - (values[0]/ values[2]));
   // printf("L1 miss rate        : %.4lf\n", values[2] * 1.0 / values[3]);
   // printf("L2 miss rate        : %.4lf\n", values[4] * 1.0 / values[2]);
    }

};




class PerfClass {
  public:

  int num_events=7;
long long *values;
  //int *Events;
  int Events[7] = {PAPI_TOT_CYC, PAPI_TOT_INS }, EventSet = PAPI_NULL;
// EventSet = PAPI_NULL;
  int retval;
  long long before_time,after_time;
    double elapsed_time;

   virtual void startMonitoring() 
   {
   // Events=(int *)calloc(num_events,sizeof(int)); 
    values=(long long *)calloc(num_events,sizeof(long long));
    retval = PAPI_create_eventset(&EventSet);
    retval = PAPI_add_events(EventSet, Events, num_events);
     before_time=PAPI_get_real_nsec();
    retval = PAPI_start(EventSet);
   }

    virtual void stopMonitoring( std::map<std::string, double> &mmap) 
    {

    retval = PAPI_stop(EventSet, values);
    after_time=PAPI_get_real_nsec();
    elapsed_time=((double)(after_time-before_time))/1.0e9;
    printf("time %.4f  PAPI_TOT_CYC %.4llu \n", elapsed_time, values[0]);
    printf("------------------------------\n");
    //printf("Hardware Counters   : %d\n", PAPI_num_counters());
   // printf("MFlops                 : %.4lf\n", values[0]/ values[1]);
   /// printf("L2 cache miss ratio        : %.4llu\n",(values[4]/values[5]) );
    printf("MIPS        : %.4f\n", values[1] / elapsed_time);
    printf("Processor utilization        : %.4f\n", (values[0] / elapsed_time));
    }
};

class PapiClass {
  public:

   // std::unique_ptr<Papi> papi_ptr; 

    int retval,cid,rapl_cid=-1,numcmp;
    int EventSet = PAPI_NULL;
    long long *values;
    int num_events=0;
    int code;
    char event_names[MAX_RAPL_EVENTS][PAPI_MAX_STR_LEN];
    char units[MAX_RAPL_EVENTS][PAPI_MIN_STR_LEN];
    int data_type[MAX_RAPL_EVENTS];
    int r,i;
    const PAPI_component_info_t *cmpinfo = NULL;
    PAPI_event_info_t evinfo;
    long long before_time,after_time;
    double elapsed_time;


    virtual void reset() {
  
    }

    virtual double getpower()
    {
//    return rapl_ptr->pkg_total_energy();
//
    };

    virtual void startPowerMonitoring() {

	/* PAPI Initialization */
     retval = PAPI_library_init( PAPI_VER_CURRENT );
     if ( retval != PAPI_VER_CURRENT ) {
	printf("PAPI_library_init failed %d \n",retval);
     }


    numcmp = PAPI_num_components();
    for(cid=0; cid<numcmp; cid++)
    {
        if ( (cmpinfo = PAPI_get_component_info(cid)) == NULL) {
        printf("PAPI_get_component_info failed\n");
        }
        if (strstr(cmpinfo->name,"rapl"))  
        {
            rapl_cid=cid;
            //printf("Found rapl component at cid %d\n",rapl_cid);
            if (cmpinfo->disabled)
            {
                printf("RAPL component disabled: %s\n", cmpinfo->disabled_reason);
            }
            break;
        }
    }
       
    /* Component not found */
     if (cid==numcmp) {
      printf("No rapl component found\n");
     }

     /* Create EventSet */
     retval = PAPI_create_eventset( &EventSet );
     /*if (retval != PAPI_OK) {
     printf("PAPI_create_eventset() %d\n",retval);
     }*/

     /* Add all events */
     code = PAPI_NATIVE_MASK;
    
     r = PAPI_enum_cmp_event( &code, PAPI_ENUM_FIRST, rapl_cid );
    //printf("beforeNum events %d %d rapl cid %d:\n", num_events,rapl_cid );
     while ( r == PAPI_OK ) {

        retval = PAPI_event_code_to_name( code, event_names[num_events] );
	if ( retval != PAPI_OK ) {
	   printf("Error translating %#x\n",code);
	   printf("PAPI_event_code_to_name %d", retval );
	}

	retval = PAPI_get_event_info(code,&evinfo);
	if (retval != PAPI_OK) {
	   printf("Error getting event info %d \n",retval);
	}

	strncpy(units[num_events],evinfo.units,sizeof(units[0])-1);
	// buffer must be null terminated to safely use strstr operation on it below
	units[num_events][sizeof(units[0])-1] = '\0';

	data_type[num_events] = evinfo.data_type;

        retval = PAPI_add_event( EventSet, code );
        if (retval != PAPI_OK) {
	  break; /* We've hit an event limit */
	}
	num_events++;

        r = PAPI_enum_cmp_event( &code, PAPI_ENUM_EVENTS, rapl_cid );
     }
//printf("afterNum events %d :\n", num_events);
     values=(long long *)calloc(num_events,sizeof(long long));
     if (values==NULL) {
	 printf("No memory %d\n",retval);
     }

     
 /* Start Counting */
     before_time=PAPI_get_real_nsec();
     retval = PAPI_start( EventSet);
     if (retval != PAPI_OK) {
	 printf("PAPI_start %d %llu \n",retval, before_time);
     }
  }

  virtual void stopPowerMonitoring( std::map<std::string, double> &mmap) 
  {
   after_time=PAPI_get_real_nsec();
    retval = PAPI_stop( EventSet, values);
    if (retval != PAPI_OK) {
	 printf("PAPI_stop() %d\n",retval);
    }

     elapsed_time=((double)(after_time-before_time))/1.0e9;

     //if (!TESTS_QUIET) {
     //   printf("\nStopping measurements, took %.3fs, gathering results...\n\n",
	 //      elapsed_time);

	//	printf("Scaled energy measurements:\n");
	//printf("Num events %d :\n", num_events);
		for(i=0;i<num_events;i++) {
		   if (strstr(units[i],"nJ")) {
			//printf("%-40s%12.6f J\t(Average Power %.1fW)\n",event_names[i],(double)values[i]/1.0e9,
			//	((double)values[i]/1.0e9)/elapsed_time);
                if(strstr(event_names[i],"rapl:::PACKAGE_ENERGY:PACKAGE"))
                {
                 mmap["PACKAGE_ENERGY"] =((double)values[i]/1.0e9)/elapsed_time;
                
                }

                 if(strstr(event_names[i],"rapl:::DRAM_ENERGY:PACKAGE"))
                {
                 mmap["DRAM_ENERGY"] =((double)values[i]/1.0e9)/elapsed_time;
                
                }

               /*  if(strstr(event_names[i],"rapl:::PSYS_ENERGY:PACKAGE"))
                {
                 mmap["PSYS_ENERGY"] =((double)values[i]/1.0e9)/elapsed_time;
               
                }*/
		   }
		}

		/*printf("\n");
		printf("Energy measurement counts:\n");

		for(i=0;i<num_events;i++) {
		   if (strstr(event_names[i],"ENERGY_CNT")) {
			  printf("%-40s%12lld\t%#08llx\n", event_names[i], values[i], values[i]);
		   }
		}*/
     //}
     
     }
  
};
 

class RaplMeasurements {
  public:
  // Reset the class so it is ready to start measuring
  std::chrono::nanoseconds m_interval;
  std::promise<void> *exitSignal = nullptr;
  std::thread *measurementThread = nullptr;

  std::unique_ptr<Rapl> rapl_ptr; 

  RaplMeasurements(std::chrono::nanoseconds interval) : m_interval(interval) {
  
  }

  virtual void reset() {
  
  }

  virtual void startPowerMonitoring() {
     /*  exitSignal = new std::promise<void>();
       std::future<void> futureObj = exitSignal->get_future();
       measurementThread = new std::thread(&RaplMeasurements::PowerMeasurementThreadFunction, this, std::move(futureObj));       
       rapl_ptr.reset( new Rapl() );   */
       
  }

  virtual void stopPowerMonitoring() {
     /*  exitSignal->set_value();
       measurementThread->join();
       delete exitSignal;   
       delete measurementThread;  */  
  };
  

 virtual double getpower()
 {
return rapl_ptr->pkg_total_energy();

 };


  virtual void PowerMeasurementThreadFunction(std::future<void> futureObj) {
     
      //std::cout << "Thread Start" << std::endl;
      //Rapl *rapl = new Rapl();
      
	    /*               // int ms_pause = 100; 
                        std::ofstream outfile ("rapl.csv", std::ios::out | std::ios::trunc);*/
	    while (futureObj.wait_for(m_interval) == std::future_status::timeout)
	    {
              rapl_ptr->sample();
                             /*   outfile << rapl->pkg_current_power() << ","
                                << rapl->pp0_current_power() << ","
                                << rapl->pp1_current_power() << ","
                                << rapl->dram_current_power() << ","
                                << rapl->total_time() << std::endl;*/
		    //  std::cout << " \t THREAD : Doing Some Work" << std::endl;
          //rapl->sample();
      }
	  // std::cout << "Thread End" << std::endl;
        std::cout << std::endl                                
                               // << "\tTotal Energy:\t" << rapl->pkg_total_energy() << " J" << std::endl
                               << "\tTotal Energy:\t" << rapl_ptr->pkg_total_energy() << " J" << std::endl
                               // << "\tTotal Dram Energy:\t" << rapl_ptr->dram_total_energy() << " J" << std::endl
                              //  << "\tAverage Power:\t" << rapl_ptr->pkg_average_power() << " W" << std::endl
                                << "\tTime:\t" << rapl_ptr->total_time() << " sec" << std::endl;
       

        //delete rapl;
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
    CacheClass cachemeasure;
    RaplMeasurements raplMeasurement;    
    PapiClass papimeasure;
    MemClass memorymeasure;
    PerfClass perfmeasure;
    

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
