struct PolHelper
{
    const double m_process_e_pol[4] = {-1., -1., 1., 1.};
    const double m_process_p_pol[4] = {-1., 1., -1., 1.};

    double e_pol(int pol_index) const
    {
        return m_process_e_pol[pol_index];
    }

    double p_pol(int pol_index) const
    {
        return m_process_p_pol[pol_index];
    }

};

PolHelper gPolHelper;
