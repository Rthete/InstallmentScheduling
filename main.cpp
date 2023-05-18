/*
 * @FilePath: \InstallmentScheduling\main.cpp
 * @Description:  
 * @Author: rthete
 * @Date: 2023-03-14 16:19:26
 * @LastEditTime: 2023-05-18 15:38:05
 */

#include "include/method.h"

int main() {
    // run_myAPMISRR(0.6, 8);
    // run_MISRR(0.3, 8); // 28152
    // run_PMIS();
    // run_APMISRR(0.6);
    // run_APMISRR_cost(0.6, 8);
    
    // test_MISRR_theta();
    // test_MISRR_all();
    // test_APMISRR_totalTime();
    // test_APMISRR_alpha();
    // test_APMISRR_beta();
    // test_APMISRR_totalTime_cost();
    // test_APMISRR_installment();
    test_myAPMISRR_installment();
    // test_MISRRL_lambda();
    // test_MISRRL_installment();
    // test_MISRRL_all();

    // compare_MISRR_and_MISRRL();
    
    return 0;
}