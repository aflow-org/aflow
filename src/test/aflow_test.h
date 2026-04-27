
#ifndef AFLOW_AFLOW_TEST_H
#define AFLOW_AFLOW_TEST_H

#include <ostream>
#include <string>

#include <AUROSTD/aurostd_xvector.h>

void PERFORM_TESTJ(std::ostream& oss);
void PERFORM_PRX(std::ostream& oss);
bool isPGM(std::string element);
void PERFORM_TEST_ALLOYS(std::ostream& oss);
void PERFORM_TEST3(std::ostream& oss);
void pinku_funcs(float x, aurostd::xvector<float>& afunc);
int pinku_main();

#endif // AFLOW_AFLOW_TEST_H
