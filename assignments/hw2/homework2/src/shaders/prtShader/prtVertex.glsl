attribute vec3 aVertexPosition;
attribute vec3 aNormalPosition;
attribute mat3 aPrecomputeLT;

uniform mat4 uModelMatrix;
uniform mat4 uViewMatrix;
uniform mat4 uProjectionMatrix;
uniform mat3 uPrecomputeLR;
uniform mat3 uPrecomputeLG;
uniform mat3 uPrecomputeLB;

varying highp vec3 vFragPos;
varying highp vec3 vNormal;
varying highp vec3 vColor;

highp float cwiseProdSum(mat3 m1, mat3 m2) {
    highp float ret = 0.0;
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            ret += m1[i][j] * m2[i][j];
        }
    }
    return ret;
}

void main(void) {

    vFragPos = (uModelMatrix * vec4(aVertexPosition, 1.0)).xyz;
    vNormal = (uModelMatrix * vec4(aNormalPosition, 0.0)).xyz;

    gl_Position = uProjectionMatrix * uViewMatrix * uModelMatrix *
                  vec4(aVertexPosition, 1.0);

    vColor = vec3(cwiseProdSum(uPrecomputeLR, aPrecomputeLT),
                  cwiseProdSum(uPrecomputeLG, aPrecomputeLT),
                  cwiseProdSum(uPrecomputeLB, aPrecomputeLT)
                  );
}