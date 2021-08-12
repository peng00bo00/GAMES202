#include "denoiser.h"

Denoiser::Denoiser() : m_useTemportal(false) {}

void Denoiser::Reprojection(const FrameInfo &frameInfo) {
    int height = m_accColor.m_height;
    int width = m_accColor.m_width;
    Matrix4x4 preWorldToScreen =
        m_preFrameInfo.m_matrix[m_preFrameInfo.m_matrix.size() - 1];
    Matrix4x4 preWorldToCamera =
        m_preFrameInfo.m_matrix[m_preFrameInfo.m_matrix.size() - 2];
    #pragma omp parallel for
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            // TODO: Reproject
            m_valid(x, y) = false;
            m_misc(x, y) = m_accColor(x, y);

            // back projection
            int curID = frameInfo.m_id(x, y);
            Float3 p = frameInfo.m_position(x, y);

            if (curID != -1) {
                Matrix4x4 M_cur_inv = Inverse(frameInfo.m_matrix[curID]);
                Matrix4x4 M_pre = m_preFrameInfo.m_matrix[curID];
                Matrix4x4 P = preWorldToScreen * M_pre * M_cur_inv;

                p = P(p, Float3::Point);

                p.x /= p.z;
                p.y /= p.z;

                // validation check
                if (p.x >= 0 && p.x < width && p.y >= 0 && p.y < height) {
                    // find ID in the previous frame
                    int i = static_cast<int>(p.x);
                    int j = static_cast<int>(p.y);
                    int preID = m_preFrameInfo.m_id(i, j);

                    if (curID == preID) {
                        m_valid(x, y) = true;
                        m_misc(x, y) = m_accColor(i, j);
                    }
                }
            }
        }
    }
    std::swap(m_misc, m_accColor);
}

void Denoiser::TemporalAccumulation(const Buffer2D<Float3> &curFilteredColor) {
    int height = m_accColor.m_height;
    int width = m_accColor.m_width;
    int kernelRadius = 3;
#pragma omp parallel for
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            // TODO: Temporal clamp
            Float3 color = m_accColor(x, y);
            Float3 mu, sig;

            // filter window
            int x1 = x - kernelRadius;
            int x2 = x + kernelRadius;
            int y1 = y - kernelRadius;
            int y2 = y + kernelRadius;

            if (x1 < 0) x1 = 0;
            if (x2 > width) x2 = width;

            if (y1 < 0) y1 = 0;
            if (y2 > height) y2 = height;

            int n = 0;
            for (size_t i = x1; i < x2; ++i)
            {
                for (size_t j = y1; j < y2; ++j)
                {
                    mu += m_accColor(i, j);
                    ++n;
                }
            }

            // mean and std
            mu = mu / n;

            for (size_t i = x1; i < x2; ++i)
            {
                for (size_t j = y1; j < y2; ++j)
                {
                    sig += Sqr(m_accColor(i, j) - mu);
                }
            }

            sig = sig / (n-1);
            sig = SafeSqrt(sig);

            color = Clamp(color, mu-sig*m_colorBoxK, mu+sig*m_colorBoxK);

            // TODO: Exponential moving average
            float alpha = 1.0;
            if (m_valid(x, y)) alpha = m_alpha;

            m_misc(x, y) = Lerp(color, curFilteredColor(x, y), alpha);
        }
    }
    std::swap(m_misc, m_accColor);
}

Buffer2D<Float3> Denoiser::Filter(const FrameInfo &frameInfo) {
    int height = frameInfo.m_beauty.m_height;
    int width = frameInfo.m_beauty.m_width;
    Buffer2D<Float3> filteredImage = CreateBuffer2D<Float3>(width, height);
    int kernelRadius = 16;
    // #pragma omp parallel for
    // for (int y = 0; y < height; y++) {
    //     for (int x = 0; x < width; x++) {
    //         // TODO: Joint bilateral filter
    //         filteredImage(x, y) = frameInfo.m_beauty(x, y);
    //         // continue;

    //         // value and weight
    //         Float3 value_sum(0.0);
    //         float w, wp, wc, wn, wd;
    //         float weight_sum= 0.0;

    //         // filter window
    //         int x1 = x - kernelRadius;
    //         int x2 = x + kernelRadius;
    //         int y1 = y - kernelRadius;
    //         int y2 = y + kernelRadius;

    //         if (x1 < 0) x1 = 0;
    //         if (x2 > width) x2 = width;

    //         if (y1 < 0) y1 = 0;
    //         if (y2 > height) y2 = height;

    //         for (size_t i = x1; i < x2; ++i)
    //         {
    //             for (size_t j = y1; j < y2; ++j)
    //             {
    //                 // weights
    //                 // position weight
    //                 wp = (i-x)*(i-x) + (j-y)*(j-y);
    //                 wp = wp / (2.0*m_sigmaCoord*m_sigmaCoord);

    //                 // color weight
    //                 Float3 Ci = frameInfo.m_beauty(x, y);
    //                 Float3 Cj = frameInfo.m_beauty(i, j);

    //                 wc = SqrDistance(Ci, Cj);
    //                 wc = wc / (2.0*m_sigmaColor*m_sigmaColor);

    //                 // normal weight
    //                 Float3 Ni = frameInfo.m_normal(x, y);
    //                 Float3 Nj = frameInfo.m_normal(i, j);

    //                 wn = SafeAcos(Dot(Ni, Nj));
    //                 wn = wn*wn / (2.0*m_sigmaNormal*m_sigmaNormal);

    //                 // plane weight
    //                 Float3 Pi = frameInfo.m_position(x, y);
    //                 Float3 Pj = frameInfo.m_position(i, j);
    //                 Float3 dP = Pj - Pi;
    //                 if (Length(dP) > 0.0) dP = Normalize(dP);

    //                 wd = Dot(Ni, dP);
    //                 wd = wd*wd / (2.0*m_sigmaPlane*m_sigmaPlane);

    //                 w = wp + wc + wn + wd;
    //                 w = exp(-w);

    //                 Float3 v = frameInfo.m_beauty(i, j);
    //                 value_sum += v * w;
    //                 weight_sum+= w;
    //             }
    //         }

    //         if (weight_sum > 0.0) filteredImage(x, y) = value_sum / weight_sum;
    //     }
    // }
    
    // À-Trous Wavelet
    int levels = 4;
    int h = 1;
    kernelRadius = 2;

    // initialize cache
    Buffer2D<Float3> cachedImage = CreateBuffer2D<Float3>(width, height);
    cachedImage.Copy(frameInfo.m_beauty);

    for (size_t l = 0; l < levels; l++)
    {
        #pragma omp parallel for
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                // value and weight
                Float3 value_sum(0.0);
                float w, wp, wc, wn, wd;
                float weight_sum= 0.0;

                // filter on current level
                for (int i = -kernelRadius; i <= kernelRadius; ++i)
                {
                    for (int j = -kernelRadius; j <= kernelRadius; ++j) {
                        int xx = x + h * i;
                        int yy = y + h * j;

                        if (xx < 0 || xx >= width || yy < 0 || yy >= height) continue;

                        // weights
                        // position weight
                        wp = (xx-x)*(xx-x) + (yy-y)*(yy-y);
                        wp = wp / (2.0*m_sigmaCoord*m_sigmaCoord);

                        // color weight
                        Float3 Ci = cachedImage(x, y);
                        Float3 Cj = cachedImage(xx, yy);

                        wc = SqrDistance(Ci, Cj);
                        wc = wc / (2.0*m_sigmaColor*m_sigmaColor);

                        // normal weight
                        Float3 Ni = frameInfo.m_normal(x, y);
                        Float3 Nj = frameInfo.m_normal(xx, yy);

                        wn = SafeAcos(Dot(Ni, Nj));
                        wn = wn*wn / (2.0*m_sigmaNormal*m_sigmaNormal);

                        // plane weight
                        Float3 Pi = frameInfo.m_position(x, y);
                        Float3 Pj = frameInfo.m_position(xx, yy);
                        Float3 dP = Pj - Pi;
                        if (Length(dP) > 0.0) dP = Normalize(dP);

                        wd = Dot(Ni, dP);
                        wd = wd*wd / (2.0*m_sigmaPlane*m_sigmaPlane);
                        
                        w = wp + wc + wn + wd;
                        w = exp(-w);

                        Float3 v = cachedImage(xx, yy);
                        value_sum += v * w;
                        weight_sum+= w;
                    }
                }

                if (weight_sum > 0.0) filteredImage(x, y) = value_sum / weight_sum;
                else filteredImage(x, y) = cachedImage(x, y);
            }
        }

        // update step length
        h *= 2;

        // update cache
        cachedImage.Copy(filteredImage);
    }
    
    return filteredImage;
}

void Denoiser::Init(const FrameInfo &frameInfo, const Buffer2D<Float3> &filteredColor) {
    m_accColor.Copy(filteredColor);
    int height = m_accColor.m_height;
    int width = m_accColor.m_width;
    m_misc = CreateBuffer2D<Float3>(width, height);
    m_valid = CreateBuffer2D<bool>(width, height);
}

void Denoiser::Maintain(const FrameInfo &frameInfo) { m_preFrameInfo = frameInfo; }

Buffer2D<Float3> Denoiser::ProcessFrame(const FrameInfo &frameInfo) {
    // Filter current frame
    Buffer2D<Float3> filteredColor;
    filteredColor = Filter(frameInfo);

    // Reproject previous frame color to current
    if (m_useTemportal) {
        Reprojection(frameInfo);
        TemporalAccumulation(filteredColor);
    } else {
        Init(frameInfo, filteredColor);
    }

    // Maintain
    Maintain(frameInfo);
    if (!m_useTemportal) {
        m_useTemportal = true;
    }
    return m_accColor;
}