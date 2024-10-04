#include "fdcl_FFTSO3.hpp"
#include "gtest/gtest.h"

namespace
{
    class FFTSO3Test : public ::testing::Test
    {
        protected:
            fdcl::FFTSO3_complex FFTSO3;
            fdcl::FFTSO3_real RFFTSO3;

            FFTSO3Test()
            {
                int l_max=5;
                FFTSO3.init(l_max);
                RFFTSO3.init(l_max);
           }

            ~FFTSO3Test() override {}
    };

    TEST_F(FFTSO3Test, ComplexTransform)
    {
        ASSERT_NEAR(FFTSO3.check_transform(), 0., 1.e-10);
    }
    TEST_F(FFTSO3Test, Weight)
    {
        ASSERT_NEAR(FFTSO3.check_weight(), 0., 1.e-10);
    }
    TEST_F(FFTSO3Test, WignerSmallD)
    {
        ASSERT_NEAR(FFTSO3.check_wigner_d(), 0., 1.e-10);
    }
    TEST_F(FFTSO3Test, DerivWignerD)
    {
        ASSERT_NEAR(FFTSO3.check_deriv_wigner_D(), 0., 1.e-10);
    }
    TEST_F(FFTSO3Test, ComplexClebschGordon)
    {
        ASSERT_NEAR(FFTSO3.check_Clebsch_Gordon(), 0., 1.e-10);
    }

    TEST_F(FFTSO3Test, WignerDReal)
    {
        ASSERT_NEAR(RFFTSO3.check_real_harmonics(), 0., 1.e-10);
    }
    TEST_F(FFTSO3Test, RealClebschGordon)
    {
        ASSERT_NEAR(RFFTSO3.check_Clebsch_Gordon(), 0., 1.e-10);
    }
    TEST_F(FFTSO3Test, DerivWignerDReal)
    {
        ASSERT_NEAR(RFFTSO3.check_deriv_real_harmonics(), 0., 1.e-10);
    }
    TEST_F(FFTSO3Test, RealTransform)
    {
        ASSERT_NEAR(RFFTSO3.check_transform(), 0., 1.e-10);
    }
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
