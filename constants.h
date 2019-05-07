#ifndef CONSTANTS_H
#define CONSTANTS_H

#define LOG2DB(x) ((x) / 2.3025851 * 20.0)
#define LOG2IN(x) ((x) / 2.3025851 * 10.0)
#define DB2LOG(x) ((x) * 2.3025851 / 20.0)
#define IN2LOG(x) ((x) * 2.3025851 / 10.0)
#define EULERGAMMA 0.57721566
#ifndef M_PI
  #define M_PI 3.14159267
#endif
#define LOGCHI2VAR (M_PI * M_PI / 6.0)
#define LOGRESBIAS 0.375

#endif
