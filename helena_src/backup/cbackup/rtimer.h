/*
 * RTimer.h
 *
 *  Created on: 05/05/2009
 *      Author: ronanrmo
 */

#ifndef RTIMER_H_
#define RTIMER_H_

#include <sys/time.h>

class RTimer{
private:
	double total_time;
	struct timeval start_time;
	struct timeval stop_time;
	int num_calls;
public:
	RTimer();
	void start();
	void zero();
	void stop();
	double get_total_time();
	int get_num_calls();
};


#endif /* RTIMER_H_ */
