/*
 * @FilePath: /InstallmentScheduling/main.cpp
 * @Description:
 * @Author: rthete
 * @Date: 2023-03-14 16:19:26
 * @LastEditTime: 2023-11-04 16:07:53
 */

#include "include/exp_1.h"
#include "include/exp_2.h"
#include "include/exp_3.h"
#include "include/method.h"

int main() {
  /**
   * @brief 调度示例测试
   *
   */
  // ModelRunner::run_MISRR(4, 4, 400, 0.2, "../data/4-servers-example/");
  // ModelRunner::run_MISRR(4, 4, 400, 0.2, "../data/4-servers-example/", {3},
  // 2);

  /**
   * @brief 容错测试
   *
   */
  // ModelRunner::run_MISRR(30, 24, 8000, 0.3, "../data/exp1-30-servers/", {6,
  // 12}, 7); ModelRunner::run_MISRR(30, 24, 8000, 0.3,
  // "../data/exp1-30-servers/", {6, 12}, 23); ModelRunner::run_MISRR(30, 24,
  // 8000, 0.3, "../data/exp1-30-servers/"); ModelRunner::run_myAPMISRR(30, 0.5,
  // 24, 8000, 0.3, "../data/exp1-30-servers/", {6, 12}, 7);

  /**
   * @brief 带容错的三个模型对比实验
   *
   */
  // exp_3::error_SIS_15();
  // exp_3::error_TolerMIS_15();
  // exp_3::error_APMISRR_15();
  // exp_3::error_SIS_30();
  // exp_3::error_TolerMIS_30();
  // exp_3::error_APMISRR_30();
  // exp_3::error_TolerMIS_30_conflict();
  // exp_3::error_cmp_3_models_ur();
  exp_3::error_save_each_time();

  // exp_2::W_vs_m_15();
  // exp_2::theta_vs_m_15();
  // exp_2::W_vs_m_30();
  // exp_2::theta_vs_m_30();

  // without_error_15_m();
  // ur_without_error_15_m();
  // without_error_15_W();
  // ur_without_error_15_W();
  // without_error_15_theta();
  // ur_without_error_15_theta();

  // without_error_30_m();
  // ur_without_error_30_m();
  // without_error_30_W();
  // ur_without_error_30_W();
  // without_error_30_theta();
  // ur_without_error_30_theta();

  // with_error_15();

  // ModelRunner::run_SIS();
  // ModelRunner::run_myAPMISRR(0.6, 8);
  // ModelRunner::run_MISRR(0.3, 8);
  // ModelRunner::run_PMIS();
  // ModelRunner::run_APMISRR(0.6);
  // ModelRunner::run_APMISRR_cost(0.6, 8);
  // ModelRunner::run_MISRRL(1, 8);
  // ModelRunner::run_MISRRLL(1, 0.7, 8);
  // ModelRunner::run_MISRRLL(0.6, 0.6, 9);

  // test_MISRR_theta();
  // test_MISRR_all();
  // test_APMISRR_totalTime();
  // test_APMISRR_alpha();
  // test_APMISRR_beta();
  // test_APMISRR_totalTime_cost();
  // test_APMISRR_installment();
  // test_myAPMISRR_installment();
  // test_MISRRL_lambda();
  // test_MISRRL_installment();
  // test_MISRRL_all();
  // test_MISRRLL_lambda2();
  // test_MISRRLL_lambda1();
  // test_MISRRLL_all();
  // test_MISRRLL_2_lambda();

  // compare_MISRR_and_MISRRL();

  return 0;
}