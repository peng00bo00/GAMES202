# GAMES202 Homework4

作者: 彭博

## 项目描述

本次作业实现了Kulla-Conty BRDF模型，主要内容包括:

- 在离线端预计算$E(\mu)$和$E_{avg}$
- 在实时端计算BRDF的能量补偿项

## 预计算$E(\mu)$

首先在离线端完成$E(\mu)$的预计算，根据文档有:

$$
\begin{aligned}
    E(\mu_0) &= \int_{0}^{2\pi} \int_{0}^{1} f_r(\mu_0, \mu_i, \phi) \mu_i d \mu_i d \phi \\
    &= \int_{0}^{2\pi} \int_{0}^{1} \frac{F(L, H) G(L, V, H) D(H)}{4 (N \cdot V) (N \cdot L) } \mu_i d \mu_i d \phi
\end{aligned}
$$

其中, $L$为入射方向，$V$为出射方向，$H$为半程向量，$N$为法向。这里取$F=1.0$，$G$和$D$的计算函数已经给出，只需完成被积函数的计算即可：

```cpp
for (int i = 0; i < sample_count; i++) {
    // TODO: To calculate (fr * ni) / p_o here
    Vec3f L = sampleList.directions[i];
    float pdf= sampleList.PDFs[i];

    float NdotL = std::max(dot(N, L), 0.0f);
    Vec3f H = normalize(L + V);

    float VoH = std::max(dot(V, H), 0.0f);

    float D = DistributionGGX(N, H, roughness);
    float G = GeometrySmith(roughness, NdotV, NdotL);

    // fr = F * D * G / (4 * NdotL * NdotV)
    // E(mu) += fr * NdotL / pdf
    float f = D * G / (4 * NdotV);

     A += f / pdf;
}
```

详细源码可参见`./lut-gen/Emu_MC.cpp`。运行`test.sh`可得到使用蒙特卡洛积分预计算的$E(\mu)$如下图所示：

<div align=center>
<img src="images/GGX_E_MC_LUT.png">
</div>

## 预计算$E_{avg}$

完成$E(\mu)$的预计算后可以在它的基础上实现$E_{avg}$的预计算。$E_{avg}$的计算公式如下：

$$
\begin{aligned}
    E_{avg} &= 2 \int_0^1 E(\mu) \mu d \mu
\end{aligned}
$$

其中$E(\mu)$可以直接读取之前的计算结果。这里同样使用蒙特卡洛积分的方式通过对立体角采样来进行计算：

```cpp
for (int i = 0; i < sample_count; i++) {
    Vec3f L = sampleList.directions[i];
    Vec3f H = normalize(V + L);

    float NoL = std::max(L.z, 0.0f);
    float NoH = std::max(H.z, 0.0f);
    float VoH = std::max(dot(V, H), 0.0f);
    float NoV = std::max(dot(N, V), 0.0f);

    // TODO: To calculate Eavg here
    float pdf= sampleList.PDFs[i];
    Eavg += Ei * NoL / (pdf * M_PI);
}
```

需要说明的是我认为这里可以进行简化，无需对立体角进行采样而是对$\mu$采样来计算$E_{avg}$。这里可以把$E_{avg}$看做是$\mu$的函数且$\mu$自身是(0, 1)范围内的变量。在主程序中已经对$\mu$进行了遍历，相当于在(0, 1)范围内进行了均匀采样，因此可以将`main()`函数修改如下：

```cpp
// Eavg += IntegrateEmu(V, roughness, NdotV, Ei) * step;
Eavg += Ei * NdotV * step * 2.0;
```

详细源码可参见`./lut-gen/Eavg_MC.cpp`。运行`test.sh`可得到使用蒙特卡洛积分预计算的$E_{avg}$如下图所示：

<div align=center>
<img src="images/GGX_Eavg_MC_LUT.png">
</div>

## Bonus1 重要性采样

## BRDF补偿

完成预计算部分后可以实现在线端的shader。首先在`PBRFragment.glsl`和`KullaContyFragment.glsl`中补充代码实现PBR材质，相关代码可以参考`./lut-gen/Emu_MC.cpp`中已经给出的实现：

```glsl
float DistributionGGX(vec3 N, vec3 H, float roughness)
{
    // TODO: To calculate GGX NDF here
    float a = roughness*roughness;
    float a2 = a*a;
    float NdotH = SaturateDot(N, H);
    float NdotH2 = NdotH*NdotH;

    float nom   = a2;
    float denom = (NdotH2 * (a2 - 1.0) + 1.0);
    denom = PI * denom * denom;

    return nom / max(denom, 0.0001);
}

float GeometrySchlickGGX(float NdotV, float roughness)
{
    // TODO: To calculate Schlick G1 here
    float a = roughness + 1.0;
    float k = (a * a) / 8.0;

    float nom = NdotV;
    float denom = NdotV * (1.0 - k) + k;

    return nom / denom;
}

float GeometrySmith(vec3 N, vec3 V, vec3 L, float roughness)
{
    // TODO: To calculate Smith G here
    float NoL = SaturateDot(N, L);
    float NoV = SaturateDot(N, V);

    float ggx2 = GeometrySchlickGGX(NoV, roughness);
    float ggx1 = GeometrySchlickGGX(NoL, roughness);

    return ggx1 * ggx2;
}

vec3 fresnelSchlick(vec3 F0, vec3 V, vec3 H)
{
    // TODO: To calculate Schlick F here
    float VoH = SaturateDot(V, H);
    return F0 + (vec3(1.0) - F0) * pow(1.0 - VoH, 5.0);
}
```

然后完成BRDF能量补偿项的计算，计算公式为：

$$
\begin{aligned}
    f_{ms} (\mu_0, \mu_i) = \frac{(1 - E(\mu_0)) (1 - E(\mu_i))}{\pi (1 - E_{avg})}
\end{aligned}
$$

$$
\begin{aligned}
    f_{add} (\mu_0, \mu_i) = \frac{F_{avg} E_{avg}}{1 - F_{avg}(1 - E_{avg})}
\end{aligned}
$$

其中$f_{ms}$为不考虑颜色情况下的能量补偿，而$f_{add}$为考虑物体自身颜色对能量吸收的情况下最终反射出能量的比例。二者相乘即为最终的能量补偿项。

在`MultiScatterBRDF()`函数中补充代码如下：
```glsl
vec3 MultiScatterBRDF(float NdotL, float NdotV)
{
    vec3 albedo = pow(texture2D(uAlbedoMap, vTextureCoord).rgb, vec3(2.2));

    vec3 E_o = texture2D(uBRDFLut, vec2(NdotL, uRoughness)).xyz;
    vec3 E_i = texture2D(uBRDFLut, vec2(NdotV, uRoughness)).xyz;

    vec3 E_avg = texture2D(uEavgLut, vec2(0, uRoughness)).xyz;
    // copper
    vec3 edgetint = vec3(0.827, 0.792, 0.678);
    vec3 F_avg = AverageFresnel(albedo, edgetint);
    
    // TODO: To calculate fms and missing energy here
    vec3 F_ms = (vec3(1.0) - E_o) * (vec3(1.0) - E_i) / (PI * (vec3(1.0) - E_avg));
    vec3 F_add= F_avg * E_avg / (vec3(1.0) - F_avg * (vec3(1.0) - E_avg));

    return F_add*F_ms;
}
```

最终得到PBR材质和Kulla-Conty材质的渲染结果如下图所示：
