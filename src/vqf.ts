const EPS = Number.EPSILON;

type Tuple<T, N extends number> = N extends N
  ? number extends N
    ? T[]
    : _TupleOf<T, N, []>
  : never;
type _TupleOf<T, N extends number, R extends unknown[]> = R["length"] extends N
  ? R
  : _TupleOf<T, N, [T, ...R]>;

export type VQFParams = {
  tauAcc: number;
  tauMag: number;
  motionBiasEstEnabled: boolean;
  restBiasEstEnabled: boolean;
  magDistRejectionEnabled: boolean;
  biasSigmaInit: number;
  biasForgettingTime: number;
  biasClip: number;
  biasSigmaMotion: number;
  biasVerticalForgettingFactor: number;
  biasSigmaRest: number;
  restMinT: number;
  restFilterTau: number;
  restThGyr: number;
  restThAcc: number;
  magCurrentTau: number;
  magRefTau: number;
  magNormTh: number;
  magDipTh: number;
  magNewTime: number;
  magNewFirstTime: number;
  magNewMinGyr: number;
  magMinUndisturbedTime: number;
  magMaxRejectionTime: number;
  magRejectionFactor: number;
};

export function defaultParams(): VQFParams {
  return {
    tauAcc: 3.0,
    tauMag: 9.0,
    motionBiasEstEnabled: true,
    restBiasEstEnabled: true,
    magDistRejectionEnabled: true,
    biasSigmaInit: 0.5,
    biasForgettingTime: 100.0,
    biasClip: 2.0,
    biasSigmaMotion: 0.1,
    biasVerticalForgettingFactor: 0.0001,
    biasSigmaRest: 0.03,
    restMinT: 1.5,
    restFilterTau: 0.5,
    restThGyr: 2.0,
    restThAcc: 0.5,
    magCurrentTau: 0.05,
    magRefTau: 20.0,
    magNormTh: 0.1,
    magDipTh: 10.0,
    magNewTime: 20.0,
    magNewFirstTime: 5.0,
    magNewMinGyr: 20.0,
    magMinUndisturbedTime: 0.5,
    magMaxRejectionTime: 60.0,
    magRejectionFactor: 2.0,
  };
}

export type VQFState = {
  gyrQuat: Tuple<number, 4>;
  accQuat: Tuple<number, 4>;
  delta: number;
  restDetected: boolean;
  magDistDetected: boolean;
  lastAccLp: Tuple<number, 3>;
  accLpState: Tuple<Tuple<number, 2>, 3>;
  lastAccCorrAngularRate: number;
  kMagInit: number;
  lastMagDisAngle: number;
  lastMagCorrAngularRate: number;
  bias: Tuple<number, 3>;
  biasP: Tuple<number, 9>;
  motionBiasEstRLpState: Tuple<Tuple<number, 2>, 9>;
  motionBiasEstBiasLpState: Tuple<Tuple<number, 2>, 2>;
  restLastSquaredDeviations: Tuple<number, 2>;
  restT: number;
  restLastGyrLp: Tuple<number, 3>;
  restGyrLpState: Tuple<Tuple<number, 2>, 3>;
  restLastAccLp: Tuple<number, 3>;
  restAccLpState: Tuple<Tuple<number, 2>, 3>;
  magRefNorm: number;
  magRefDip: number;
  magUndisturbedT: number;
  magRejectT: number;
  magCandidateNorm: number;
  magCandidateDip: number;
  magCandidateT: number;
  magNormDip: Tuple<number, 2>;
  magNormDipLpState: Tuple<Tuple<number, 2>, 2>;
};

function tupleOfSize<N extends number, T>(
  length: N,
  valueFn: (i: number) => T,
): Tuple<T, N> {
  return [...Array(length)].map((_, i) => valueFn(i)) as Tuple<T, N>;
}

function emptyQuat(): Tuple<number, 4> {
  return [0, 0, 0, 0];
}

function defaultState(): VQFState {
  return {
    gyrQuat: tupleOfSize(4, () => 0),
    accQuat: tupleOfSize(4, () => 0),
    delta: 0,
    restDetected: false,
    magDistDetected: false,
    lastAccLp: tupleOfSize(3, () => 0),
    accLpState: tupleOfSize(3, () => tupleOfSize(2, () => 0)),
    lastAccCorrAngularRate: 0,
    kMagInit: 0,
    lastMagDisAngle: 0,
    lastMagCorrAngularRate: 0,
    bias: tupleOfSize(3, () => 0),
    biasP: tupleOfSize(9, () => 0),
    motionBiasEstRLpState: tupleOfSize(9, () => tupleOfSize(2, () => 0)),
    motionBiasEstBiasLpState: tupleOfSize(2, () => tupleOfSize(2, () => 0)),
    restLastSquaredDeviations: tupleOfSize(2, () => 0),
    restT: 0,
    restLastGyrLp: tupleOfSize(3, () => 0),
    restGyrLpState: tupleOfSize(3, () => tupleOfSize(2, () => 0)),
    restLastAccLp: tupleOfSize(3, () => 0),
    restAccLpState: tupleOfSize(3, () => tupleOfSize(2, () => 0)),
    magRefNorm: 0,
    magRefDip: 0,
    magUndisturbedT: 0,
    magRejectT: 0,
    magCandidateNorm: 0,
    magCandidateDip: 0,
    magCandidateT: 0,
    magNormDip: tupleOfSize(2, () => 0),
    magNormDipLpState: tupleOfSize(2, () => tupleOfSize(2, () => 0)),
  };
}

export type VQFCoefficients = {
  gyrTs: number;
  accTs: number;
  magTs: number;
  accLpB: Tuple<number, 3>;
  accLpA: Tuple<number, 2>;
  kMag: number;
  biasP0: number;
  biasV: number;
  biasMotionW: number;
  biasVerticalW: number;
  biasRestW: number;
  restGyrLpB: Tuple<number, 3>;
  restGyrLpA: Tuple<number, 2>;
  restAccLpB: Tuple<number, 3>;
  restAccLpA: Tuple<number, 2>;
  kMagRef: number;
  magNormDipLpB: Tuple<number, 3>;
  magNormDipLpA: Tuple<number, 2>;
};

function defaultCoeffs(): VQFCoefficients {
  return {
    gyrTs: 0,
    accTs: 0,
    magTs: 0,
    accLpB: tupleOfSize(3, () => 0),
    accLpA: tupleOfSize(2, () => 0),
    kMag: 0,
    biasP0: 0,
    biasV: 0,
    biasMotionW: 0,
    biasVerticalW: 0,
    biasRestW: 0,
    restGyrLpB: tupleOfSize(3, () => 0),
    restGyrLpA: tupleOfSize(2, () => 0),
    restAccLpB: tupleOfSize(3, () => 0),
    restAccLpA: tupleOfSize(2, () => 0),
    kMagRef: 0,
    magNormDipLpB: tupleOfSize(3, () => 0),
    magNormDipLpA: tupleOfSize(2, () => 0),
  };
}

export type UpdateBatchResult<N extends number> = {
  out6D: Tuple<Tuple<number, 4>, N>;
  out9D: Tuple<Tuple<number, 4>, N>;
  outDelta: Tuple<number, N>;
  outBias: Tuple<Tuple<number, 3>, N>;
  outBiasSigma: Tuple<number, N>;
  outRest: Tuple<boolean, N>;
  outMagDist: Tuple<boolean, N>;
};

function square(x: number): number {
  return x * x;
}

export class VQF {
  constructor(
    gyrTs: number,
    accTs: number = -1,
    magTs: number = -1,
    params?: VQFParams,
  ) {
    if (params) {
      this.params = params;
    }

    this.coeffs.gyrTs = gyrTs;
    this.coeffs.accTs = accTs > 0 ? accTs : gyrTs;
    this.coeffs.magTs = magTs > 0 ? magTs : gyrTs;

    this.setup();
  }

  public updateGyr(gyr: Tuple<number, 3>): void {
    if (this.params.restBiasEstEnabled || this.params.magDistRejectionEnabled) {
      const result = VQF.filterVec(
        gyr,
        3,
        this.params.restFilterTau,
        this.coeffs.gyrTs,
        this.coeffs.restGyrLpB,
        this.coeffs.restGyrLpA,
        this.state.restGyrLpState,
      );
      this.state.restGyrLpState = result.state;
      this.state.restLastGyrLp = result.out;

      this.state.restLastSquaredDeviations[0] =
        square(gyr[0] - this.state.restLastGyrLp[0]) +
        square(gyr[1] - this.state.restLastGyrLp[1]) +
        square(gyr[2] - this.state.restLastGyrLp[2]);

      const biasClip = this.params.biasClip * (Math.PI / 180.0);
      if (
        this.state.restLastSquaredDeviations[0] >=
          square(this.params.restThGyr * (Math.PI / 180.0)) ||
        Math.abs(this.state.restLastGyrLp[0]) > biasClip ||
        Math.abs(this.state.restLastGyrLp[1]) > biasClip ||
        Math.abs(this.state.restLastGyrLp[2]) > biasClip
      ) {
        this.state.restT = 0.0;
        this.state.restDetected = false;
      }
    }

    // remove estimated gyro bias
    const gyrNoBias: Tuple<number, 3> = [
      gyr[0] - this.state.bias[0],
      gyr[1] - this.state.bias[1],
      gyr[2] - this.state.bias[2],
    ];

    // gyroscope prediction step
    const gyrNorm = VQF.norm(gyrNoBias, 3);
    const angle = gyrNorm * this.coeffs.gyrTs;
    if (gyrNorm > EPS) {
      const c = Math.cos(angle / 2);
      const s = Math.sin(angle / 2) / gyrNorm;
      const gyrStepQuat: Tuple<number, 4> = [
        c,
        s * gyrNoBias[0],
        s * gyrNoBias[1],
        s * gyrNoBias[2],
      ];
      this.state.gyrQuat = VQF.quatMultiply(this.state.gyrQuat, gyrStepQuat);
      this.state.gyrQuat = VQF.normalize(this.state.gyrQuat, 4);
    }
  }

  public updateAcc(acc: Tuple<number, 3>): void {
    if (acc[0] == 0.0 && acc[1] == 0.0 && acc[2] == 0.0) {
      return;
    }

    if (this.params.restBiasEstEnabled) {
      const result = VQF.filterVec(
        acc,
        3,
        this.params.restFilterTau,
        this.coeffs.accTs,
        this.coeffs.restAccLpB,
        this.coeffs.restAccLpA,
        this.state.restAccLpState,
      );

      this.state.restAccLpState = result.state;
      this.state.restLastAccLp = result.out;

      this.state.restLastSquaredDeviations[1] =
        square(acc[0] - this.state.restLastAccLp[0]) +
        square(acc[1] - this.state.restLastAccLp[1]) +
        square(acc[2] - this.state.restLastAccLp[2]);

      if (
        this.state.restLastSquaredDeviations[1] >= square(this.params.restThAcc)
      ) {
        this.state.restT = 0.0;
        this.state.restDetected = false;
      } else {
        this.state.restT += this.coeffs.accTs;
        if (this.state.restT >= this.params.restMinT) {
          this.state.restDetected = true;
        }
      }
    }

    let accEarth = tupleOfSize(3, () => 0);

    // filter acc in inertial frame
    accEarth = VQF.quatRotate(this.state.gyrQuat, acc);
    const result = VQF.filterVec(
      accEarth,
      3,
      this.params.tauAcc,
      this.coeffs.accTs,
      this.coeffs.accLpB,
      this.coeffs.accLpA,
      this.state.accLpState,
    );
    this.state.accLpState = result.state;
    this.state.lastAccLp = result.out;

    // transform to 6D earth frame and normalize
    accEarth = VQF.quatRotate(this.state.accQuat, this.state.lastAccLp);
    accEarth = VQF.normalize(accEarth, 3);

    // inclination correction
    const accCorrQuat: Tuple<number, 4> = [0, 0, 0, 0];
    const q_w = Math.sqrt((accEarth[2] + 1) / 2);
    if (q_w > 1e-6) {
      accCorrQuat[0] = q_w;
      accCorrQuat[1] = (0.5 * accEarth[1]) / q_w;
      accCorrQuat[2] = (-0.5 * accEarth[0]) / q_w;
      accCorrQuat[3] = 0;
    } else {
      // to avoid numeric issues when acc is close to [0 0 -1], i.e. the correction step is close (<= 0.00011°) to 180°:
      accCorrQuat[0] = 0;
      accCorrQuat[1] = 1;
      accCorrQuat[2] = 0;
      accCorrQuat[3] = 0;
    }
    this.state.accQuat = VQF.quatMultiply(accCorrQuat, this.state.accQuat);
    this.state.accQuat = VQF.normalize(this.state.accQuat, 4);

    // calculate correction angular rate to facilitate debugging
    this.state.lastAccCorrAngularRate =
      Math.acos(accEarth[2]) / this.coeffs.accTs;

    // bias estimation
    if (this.params.motionBiasEstEnabled || this.params.restBiasEstEnabled) {
      let biasClip = this.params.biasClip * (Math.PI / 180.0);

      let accGyrQuat = emptyQuat();
      let R = tupleOfSize(9, () => 0);
      let biasLp = tupleOfSize(2, () => 0);

      // get rotation matrix corresponding to accGyrQuat
      accGyrQuat = this.getQuat6D();
      R[0] = 1 - 2 * square(accGyrQuat[2]) - 2 * square(accGyrQuat[3]); // r11
      R[1] =
        2 * (accGyrQuat[2] * accGyrQuat[1] - accGyrQuat[0] * accGyrQuat[3]); // r12
      R[2] =
        2 * (accGyrQuat[0] * accGyrQuat[2] + accGyrQuat[3] * accGyrQuat[1]); // r13
      R[3] =
        2 * (accGyrQuat[0] * accGyrQuat[3] + accGyrQuat[2] * accGyrQuat[1]); // r21
      R[4] = 1 - 2 * square(accGyrQuat[1]) - 2 * square(accGyrQuat[3]); // r22
      R[5] =
        2 * (accGyrQuat[2] * accGyrQuat[3] - accGyrQuat[1] * accGyrQuat[0]); // r23
      R[6] =
        2 * (accGyrQuat[3] * accGyrQuat[1] - accGyrQuat[0] * accGyrQuat[2]); // r31
      R[7] =
        2 * (accGyrQuat[0] * accGyrQuat[1] + accGyrQuat[3] * accGyrQuat[2]); // r32
      R[8] = 1 - 2 * square(accGyrQuat[1]) - 2 * square(accGyrQuat[2]); // r33

      // calculate R*b_hat (only the x and y component, as z is not needed)
      biasLp[0] =
        R[0] * this.state.bias[0] +
        R[1] * this.state.bias[1] +
        R[2] * this.state.bias[2];
      biasLp[1] =
        R[3] * this.state.bias[0] +
        R[4] * this.state.bias[1] +
        R[5] * this.state.bias[2];

      // low-pass filter R and R*b_hat
      const filteredR = VQF.filterVec(
        R,
        9,
        this.params.tauAcc,
        this.coeffs.accTs,
        this.coeffs.accLpB,
        this.coeffs.accLpA,
        this.state.motionBiasEstRLpState,
      );
      this.state.motionBiasEstRLpState = filteredR.state;
      R = filteredR.out;

      const filteredRbhat = VQF.filterVec(
        biasLp,
        2,
        this.params.tauAcc,
        this.coeffs.accTs,
        this.coeffs.accLpB,
        this.coeffs.accLpA,
        this.state.motionBiasEstBiasLpState,
      );
      this.state.motionBiasEstBiasLpState = filteredRbhat.state;
      biasLp = filteredRbhat.out;

      // set measurement error and covariance for the respective Kalman filter update
      let w = tupleOfSize(3, () => 0);
      let e = tupleOfSize(3, () => 0);
      if (this.state.restDetected && this.params.restBiasEstEnabled) {
        e[0] = this.state.restLastGyrLp[0] - this.state.bias[0];
        e[1] = this.state.restLastGyrLp[1] - this.state.bias[1];
        e[2] = this.state.restLastGyrLp[2] - this.state.bias[2];
        R = VQF.matrix3SetToScaledIdentity(1.0);
        w.fill(this.coeffs.biasRestW);
      } else if (this.params.motionBiasEstEnabled) {
        e[0] =
          -accEarth[1] / this.coeffs.accTs +
          biasLp[0] -
          R[0] * this.state.bias[0] -
          R[1] * this.state.bias[1] -
          R[2] * this.state.bias[2];
        e[1] =
          accEarth[0] / this.coeffs.accTs +
          biasLp[1] -
          R[3] * this.state.bias[0] -
          R[4] * this.state.bias[1] -
          R[5] * this.state.bias[2];
        e[2] =
          -R[6] * this.state.bias[0] -
          R[7] * this.state.bias[1] -
          R[8] * this.state.bias[2];
        w[0] = this.coeffs.biasMotionW;
        w[1] = this.coeffs.biasMotionW;
        w[2] = this.coeffs.biasVerticalW;
      } else {
        w.fill(-1);
      }

      // Kalman filter update
      // step 1: P = P + V (also increase covariance if there is no measurement update!)
      if (this.state.biasP[0] < this.coeffs.biasP0) {
        this.state.biasP[0] += this.coeffs.biasV;
      }
      if (this.state.biasP[4] < this.coeffs.biasP0) {
        this.state.biasP[4] += this.coeffs.biasV;
      }
      if (this.state.biasP[8] < this.coeffs.biasP0) {
        this.state.biasP[8] += this.coeffs.biasV;
      }
      if (w[0] >= 0) {
        // clip disagreement to -2..2 °/s
        // (this also effectively limits the harm done by the first inclination correction step)
        e = VQF.clip(e, -biasClip, biasClip, 3);

        // step 2: K = P R^T inv(W + R P R^T)
        let K = tupleOfSize(9, () => 0);
        K = VQF.matrix3MultiplyTpsSecond(this.state.biasP, R); // K = P R^T
        K = VQF.matrix3Multiply(R, K); // K = R P R^T
        K[0] += w[0];
        K[4] += w[1];
        K[8] += w[2]; // K = W + R P R^T
        K = VQF.matrix3Inv(K); // K = inv(W + R P R^T)
        K = VQF.matrix3MultiplyTpsFirst(R, K); // K = R^T inv(W + R P R^T)
        K = VQF.matrix3Multiply(this.state.biasP, K); // K = P R^T inv(W + R P R^T)

        // step 3: bias = bias + K (y - R bias) = bias + K e
        this.state.bias[0] += K[0] * e[0] + K[1] * e[1] + K[2] * e[2];
        this.state.bias[1] += K[3] * e[0] + K[4] * e[1] + K[5] * e[2];
        this.state.bias[2] += K[6] * e[0] + K[7] * e[1] + K[8] * e[2];

        // step 4: P = P - K R P
        K = VQF.matrix3Multiply(K, R); // K = K R
        K = VQF.matrix3Multiply(K, this.state.biasP); // K = K R P
        for (let i = 0; i < 9; i++) {
          this.state.biasP[i] -= K[i];
        }

        // clip bias estimate to -2..2 °/s
        this.state.bias = VQF.clip(this.state.bias, -biasClip, biasClip, 3);
      }
    }
  }

  public updateMag(mag: Tuple<number, 3>): void {
    // ignore [0 0 0] samples
    if (mag[0] == 0.0 && mag[1] == 0.0 && mag[2] == 0.0) {
      return;
    }

    let magEarth = tupleOfSize(3, () => 0);

    // bring magnetometer measurement into 6D earth frame
    let accGyrQuat = emptyQuat();
    accGyrQuat = this.getQuat6D();
    magEarth = VQF.quatRotate(accGyrQuat, mag);

    if (this.params.magDistRejectionEnabled) {
      this.state.magNormDip[0] = VQF.norm(magEarth, 3);
      this.state.magNormDip[1] = -Math.asin(
        magEarth[2] / this.state.magNormDip[0],
      );

      if (this.params.magCurrentTau > 0) {
        const result = VQF.filterVec(
          this.state.magNormDip,
          2,
          this.params.magCurrentTau,
          this.coeffs.magTs,
          this.coeffs.magNormDipLpB,
          this.coeffs.magNormDipLpA,
          this.state.magNormDipLpState,
        );
        this.state.magNormDipLpState = result.state;
        this.state.magNormDip = result.out;
      }

      // magnetic disturbance detection
      if (
        Math.abs(this.state.magNormDip[0] - this.state.magRefNorm) <
          this.params.magNormTh * this.state.magRefNorm &&
        Math.abs(this.state.magNormDip[1] - this.state.magRefDip) <
          this.params.magDipTh * (Math.PI / 180.0)
      ) {
        this.state.magUndisturbedT += this.coeffs.magTs;
        if (this.state.magUndisturbedT >= this.params.magMinUndisturbedTime) {
          this.state.magDistDetected = false;
          this.state.magRefNorm +=
            this.coeffs.kMagRef *
            (this.state.magNormDip[0] - this.state.magRefNorm);
          this.state.magRefDip +=
            this.coeffs.kMagRef *
            (this.state.magNormDip[1] - this.state.magRefDip);
        }
      } else {
        this.state.magUndisturbedT = 0.0;
        this.state.magDistDetected = true;
      }

      // new magnetic field acceptance
      if (
        Math.abs(this.state.magNormDip[0] - this.state.magCandidateNorm) <
          this.params.magNormTh * this.state.magCandidateNorm &&
        Math.abs(this.state.magNormDip[1] - this.state.magCandidateDip) <
          this.params.magDipTh * (Math.PI / 180.0)
      ) {
        if (
          VQF.norm(this.state.restLastGyrLp, 3) >=
          (this.params.magNewMinGyr * Math.PI) / 180.0
        ) {
          this.state.magCandidateT += this.coeffs.magTs;
        }
        this.state.magCandidateNorm +=
          this.coeffs.kMagRef *
          (this.state.magNormDip[0] - this.state.magCandidateNorm);
        this.state.magCandidateDip +=
          this.coeffs.kMagRef *
          (this.state.magNormDip[1] - this.state.magCandidateDip);

        if (
          this.state.magDistDetected &&
          (this.state.magCandidateT >= this.params.magNewTime ||
            (this.state.magRefNorm == 0.0 &&
              this.state.magCandidateT >= this.params.magNewFirstTime))
        ) {
          this.state.magRefNorm = this.state.magCandidateNorm;
          this.state.magRefDip = this.state.magCandidateDip;
          this.state.magDistDetected = false;
          this.state.magUndisturbedT = this.params.magMinUndisturbedTime;
        }
      } else {
        this.state.magCandidateT = 0.0;
        this.state.magCandidateNorm = this.state.magNormDip[0];
        this.state.magCandidateDip = this.state.magNormDip[1];
      }
    }

    // calculate disagreement angle based on current magnetometer measurement
    this.state.lastMagDisAngle =
      Math.atan2(magEarth[0], magEarth[1]) - this.state.delta;

    // make sure the disagreement angle is in the range [-pi, pi]
    if (this.state.lastMagDisAngle > Math.PI) {
      this.state.lastMagDisAngle -= 2 * Math.PI;
    } else if (this.state.lastMagDisAngle < -Math.PI) {
      this.state.lastMagDisAngle += 2 * Math.PI;
    }

    let k = this.coeffs.kMag;

    if (this.params.magDistRejectionEnabled) {
      // magnetic disturbance rejection
      if (this.state.magDistDetected) {
        if (this.state.magRejectT <= this.params.magMaxRejectionTime) {
          this.state.magRejectT += this.coeffs.magTs;
          k = 0;
        } else {
          k /= this.params.magRejectionFactor;
        }
      } else {
        this.state.magRejectT = Math.max(
          this.state.magRejectT -
            this.params.magRejectionFactor * this.coeffs.magTs,
          0.0,
        );
      }
    }

    // ensure fast initial convergence
    if (this.state.kMagInit != 0.0) {
      // make sure that the gain k is at least 1/N, N=1,2,3,... in the first few samples
      if (k < this.state.kMagInit) {
        k = this.state.kMagInit;
      }

      // iterative expression to calculate 1/N
      this.state.kMagInit = this.state.kMagInit / (this.state.kMagInit + 1);

      // disable if t > tauMag
      if (this.state.kMagInit * this.params.tauMag < this.coeffs.magTs) {
        this.state.kMagInit = 0.0;
      }
    }

    // first-order filter step
    this.state.delta += k * this.state.lastMagDisAngle;
    // calculate correction angular rate to facilitate debugging
    this.state.lastMagCorrAngularRate =
      (k * this.state.lastMagDisAngle) / this.coeffs.magTs;

    // make sure delta is in the range [-pi, pi]
    if (this.state.delta > Math.PI) {
      this.state.delta -= 2 * Math.PI;
    } else if (this.state.delta < -Math.PI) {
      this.state.delta += 2 * Math.PI;
    }
  }

  public update(
    gyr: Tuple<number, 3>,
    acc: Tuple<number, 3>,
    mag?: Tuple<number, 3>,
  ): void {
    this.updateGyr(gyr);
    this.updateAcc(acc);
    if (mag) {
      this.updateMag(mag);
    }
  }

  public updateBatch<N extends number>(
    gyr: Tuple<Tuple<number, 3>, N>,
    acc: Tuple<Tuple<number, 3>, N>,
    mag: Tuple<Tuple<number, 3>, N> | undefined,
    n: N,
  ): UpdateBatchResult<N> {
    const out6D = tupleOfSize(length, () => emptyQuat());
    const out9D = tupleOfSize(length, () => emptyQuat());
    const outDelta = tupleOfSize(length, () => 0);
    const outBias = tupleOfSize(length, () => tupleOfSize(3, () => 0));
    const outBiasSigma = tupleOfSize(length, () => 0);
    const outRest = tupleOfSize(length, () => false);
    const outMagDist = tupleOfSize(length, () => false);

    for (let i = 0; i < n; i++) {
      if (mag) {
        this.update(gyr[i], acc[i], mag[i]);
      } else {
        this.update(gyr[i], acc[i]);
      }
      out6D[i] = this.getQuat6D();
      out9D[i] = this.getQuat9D();
      outDelta[i] = this.state.delta;
      outBias[i] = this.state.bias;
      outBiasSigma[i] = this.getBiasEstimate()[0];
      outRest[i] = this.state.restDetected;
      outMagDist[i] = this.state.magDistDetected;
    }

    return {
      out6D: out6D as Tuple<Tuple<number, 4>, N>,
      out9D: out9D as Tuple<Tuple<number, 4>, N>,
      outDelta: outDelta as Tuple<number, N>,
      outBias: outBias as Tuple<Tuple<number, 3>, N>,
      outBiasSigma: outBiasSigma as Tuple<number, N>,
      outRest: outRest as Tuple<boolean, N>,
      outMagDist: outMagDist as Tuple<boolean, N>,
    };
  }

  public getQuat3D(): Tuple<number, 4> {
    return this.state.gyrQuat;
  }

  public getQuat6D(): Tuple<number, 4> {
    return VQF.quatMultiply(this.state.accQuat, this.state.gyrQuat);
  }

  public getQuat9D(): Tuple<number, 4> {
    return VQF.quatApplyDelta(this.getQuat6D(), this.state.delta);
  }

  public getDelta(): number {
    return this.state.delta;
  }

  public getBiasEstimate(): [number, Tuple<number, 3>] {
    // use largest absolute row sum as upper bound estimate for largest eigenvalue (Gershgorin circle theorem)
    // and clip output to biasSigmaInit
    const sum1 =
      Math.abs(this.state.biasP[0]) +
      Math.abs(this.state.biasP[1]) +
      Math.abs(this.state.biasP[2]);
    const sum2 =
      Math.abs(this.state.biasP[3]) +
      Math.abs(this.state.biasP[4]) +
      Math.abs(this.state.biasP[5]);
    const sum3 =
      Math.abs(this.state.biasP[6]) +
      Math.abs(this.state.biasP[7]) +
      Math.abs(this.state.biasP[8]);
    const P = Math.min(Math.max(sum1, sum2, sum3), this.coeffs.biasP0);

    return [Math.sqrt(P) * (Math.PI / 100 / 180), this.state.bias];
  }

  public setBiasEstimate(bias: Tuple<number, 3>, sigma: number) {
    this.state.bias = bias;
    if (sigma > 0) {
      const P = square(sigma * ((180.0 * 100.0) / Math.PI));
      this.state.biasP = VQF.matrix3SetToScaledIdentity(P);
    }
  }

  public getRestDetected(): boolean {
    return this.state.restDetected;
  }

  public getMagDistDetected(): boolean {
    return this.state.magDistDetected;
  }

  public getRelativeRestDeviations(): Tuple<number, 2> {
    return [
      Math.sqrt(this.state.restLastSquaredDeviations[0]) /
        (this.params.restThGyr * (Math.PI / 180.0)),
      Math.sqrt(this.state.restLastSquaredDeviations[1]) /
        this.params.restThAcc,
    ];
  }

  public getMagRefNorm(): number {
    return this.state.magRefNorm;
  }

  public getMagRefDip(): number {
    return this.state.magRefDip;
  }

  public setMagRef(norm: number, dip: number): void {
    this.state.magRefNorm = norm;
    this.state.magRefDip = dip;
  }

  public setTauAcc(tauAcc: number): void {
    if (this.params.tauAcc == tauAcc) {
      return;
    }
    this.params.tauAcc = tauAcc;

    let [newB, newA] = VQF.filterCoeffs(this.params.tauAcc, this.coeffs.accTs);
    this.state.accLpState = VQF.filterAdaptStateForCoeffChange(
      this.state.lastAccLp,
      3,
      this.coeffs.accLpB,
      this.coeffs.accLpA,
      newB,
      newA,
      this.state.accLpState,
    );

    // For R and biasLP, the last value is not saved in the state.
    // Since b0 is small (at reasonable settings), the last output is close to state[0].
    const R = tupleOfSize(9, () => 0);
    for (let i = 0; i < 9; i++) {
      R[i] = this.state.motionBiasEstRLpState[i][0];
    }
    this.state.motionBiasEstRLpState = VQF.filterAdaptStateForCoeffChange(
      R,
      9,
      this.coeffs.accLpB,
      this.coeffs.accLpA,
      newB,
      newA,
      this.state.motionBiasEstRLpState,
    );
    const biasLp = tupleOfSize(2, () => 0);
    for (let i = 0; i < 2; i++) {
      biasLp[i] = this.state.motionBiasEstBiasLpState[i][0];
    }
    this.state.motionBiasEstBiasLpState = VQF.filterAdaptStateForCoeffChange(
      biasLp,
      2,
      this.coeffs.accLpB,
      this.coeffs.accLpA,
      newB,
      newA,
      this.state.motionBiasEstBiasLpState,
    );

    this.coeffs.accLpB = newB;
    this.coeffs.accLpA = newA;
  }

  public setTauMag(tauMag: number): void {
    this.params.tauMag = tauMag;
    this.coeffs.kMag = VQF.gainFromTau(this.params.tauMag, this.coeffs.magTs);
  }

  public setMotionBiasEstEnabled(enabled: boolean): void {
    if (this.params.motionBiasEstEnabled == enabled) {
      return;
    }
    this.params.motionBiasEstEnabled = enabled;
    this.state.motionBiasEstRLpState = tupleOfSize(9, () =>
      tupleOfSize(2, () => NaN),
    );
    this.state.motionBiasEstBiasLpState = tupleOfSize(2, () =>
      tupleOfSize(2, () => NaN),
    );
  }

  public setRestBiasEstEnabled(enabled: boolean): void {
    if (this.params.restBiasEstEnabled == enabled) {
      return;
    }
    this.params.restBiasEstEnabled = enabled;
    this.state.restDetected = false;
    this.state.restLastSquaredDeviations.fill(0);
    this.state.restT = 0.0;
    this.state.restLastGyrLp.fill(0);
    this.state.restGyrLpState = tupleOfSize(3, () => tupleOfSize(2, () => NaN));
    this.state.restLastAccLp.fill(0);
    this.state.restAccLpState = tupleOfSize(3, () => tupleOfSize(2, () => NaN));
  }

  public setMagDistRejectionEnabled(enabled: boolean): void {
    if (this.params.magDistRejectionEnabled == enabled) {
      return;
    }
    this.params.magDistRejectionEnabled = enabled;
    this.state.magDistDetected = true;
    this.state.magRefNorm = 0.0;
    this.state.magRefDip = 0.0;
    this.state.magUndisturbedT = 0.0;
    this.state.magRejectT = this.params.magMaxRejectionTime;
    this.state.magCandidateNorm = -1.0;
    this.state.magCandidateDip = 0.0;
    this.state.magCandidateT = 0.0;
    this.state.magNormDipLpState = tupleOfSize(2, () =>
      tupleOfSize(2, () => NaN),
    );
  }

  public setRestDetectionThresholds(thGyr: number, thAcc: number): void {
    this.params.restThGyr = thGyr;
    this.params.restThAcc = thAcc;
  }

  public getParams(): VQFParams {
    return this.params;
  }

  public getCoeffs(): VQFCoefficients {
    return this.coeffs;
  }

  public getState(): VQFState {
    return this.state;
  }

  public setState(state: VQFState): void {
    this.state = state;
  }

  public resetState(): void {
    this.state.gyrQuat = VQF.quatSetToIdentity();
    this.state.accQuat = VQF.quatSetToIdentity();
    this.state.delta = 0.0;

    this.state.restDetected = false;
    this.state.magDistDetected = true;

    this.state.lastAccLp.fill(0);
    this.state.accLpState = tupleOfSize(3, () => tupleOfSize(2, () => 0));
    this.state.lastAccCorrAngularRate = 0.0;

    this.state.kMagInit = 1.0;
    this.state.lastMagDisAngle = 0.0;
    this.state.lastMagCorrAngularRate = 0.0;

    this.state.bias.fill(0);
    this.state.biasP = VQF.matrix3SetToScaledIdentity(this.coeffs.biasP0);

    this.state.motionBiasEstRLpState = tupleOfSize(9, () =>
      tupleOfSize(2, () => NaN),
    );
    this.state.motionBiasEstBiasLpState = tupleOfSize(2, () =>
      tupleOfSize(2, () => NaN),
    );

    this.state.restLastSquaredDeviations.fill(0);
    this.state.restT = 0.0;
    this.state.restLastGyrLp.fill(0);
    this.state.restGyrLpState = tupleOfSize(3, () => tupleOfSize(2, () => NaN));
    this.state.restLastAccLp.fill(0);
    this.state.restAccLpState = tupleOfSize(3, () => tupleOfSize(2, () => NaN));

    this.state.magRefNorm = 0.0;
    this.state.magRefDip = 0.0;
    this.state.magUndisturbedT = 0.0;
    this.state.magRejectT = this.params.magMaxRejectionTime;
    this.state.magCandidateNorm = -1.0;
    this.state.magCandidateDip = 0.0;
    this.state.magCandidateT = 0.0;
    this.state.magNormDip.fill(0);
    this.state.magNormDipLpState = tupleOfSize(2, () =>
      tupleOfSize(2, () => NaN),
    );
  }

  public static quatMultiply(
    q1: Tuple<number, 4>,
    q2: Tuple<number, 4>,
  ): Tuple<number, 4> {
    const w = q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3];
    const x = q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2];
    const y = q1[0] * q2[2] - q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1];
    const z = q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0];
    return [w, x, y, z];
  }

  public static quatConj(q: Tuple<number, 4>): Tuple<number, 4> {
    const w = q[0];
    const x = -q[1];
    const y = -q[2];
    const z = -q[3];
    return [w, x, y, z];
  }

  public static quatSetToIdentity(): Tuple<number, 4> {
    return [1, 0, 0, 0];
  }

  public static quatApplyDelta(
    q: Tuple<number, 4>,
    delta: number,
  ): Tuple<number, 4> {
    // out = quatMultiply([cos(delta/2), 0, 0, sin(delta/2)], q)
    const c = Math.cos(delta / 2);
    const s = Math.sin(delta / 2);
    const w = c * q[0] - s * q[3];
    const x = c * q[1] - s * q[2];
    const y = c * q[2] + s * q[1];
    const z = c * q[3] + s * q[0];
    return [x, y, z, w];
  }

  public static quatRotate(
    q: Tuple<number, 4>,
    v: Tuple<number, 3>,
  ): Tuple<number, 3> {
    const x =
      (1 - 2 * q[2] * q[2] - 2 * q[3] * q[3]) * v[0] +
      2 * v[1] * (q[2] * q[1] - q[0] * q[3]) +
      2 * v[2] * (q[0] * q[2] + q[3] * q[1]);
    const y =
      2 * v[0] * (q[0] * q[3] + q[2] * q[1]) +
      v[1] * (1 - 2 * q[1] * q[1] - 2 * q[3] * q[3]) +
      2 * v[2] * (q[2] * q[3] - q[1] * q[0]);
    const z =
      2 * v[0] * (q[3] * q[1] - q[0] * q[2]) +
      2 * v[1] * (q[0] * q[1] + q[3] * q[2]) +
      v[2] * (1 - 2 * q[1] * q[1] - 2 * q[2] * q[2]);
    return [x, y, z];
  }

  public static norm<N extends number>(vec: Tuple<number, N>, n: N): number {
    let s = 0;
    for (let i = 0; i < n; i++) {
      s += vec[i] * vec[i];
    }
    return Math.sqrt(s);
  }

  public static normalize<N extends number>(
    vec: Tuple<number, N>,
    N: N,
  ): Tuple<number, N> {
    const n = VQF.norm(vec, N);
    if (n < EPS) {
      return vec;
    }
    for (let i = 0; i < N; i++) {
      vec[i] /= n;
    }

    return vec;
  }

  public static clip<N extends number>(
    vec: Tuple<number, N>,
    min: number,
    max: number,
    N: N,
  ): Tuple<number, N> {
    for (let i = 0; i < N; i++) {
      if (vec[i] < min) {
        vec[i] = min;
      } else if (vec[i] > max) {
        vec[i] = max;
      }
    }
    return vec;
  }

  public static gainFromTau(tau: number, Ts: number): number {
    if (tau < 0.0) {
      return 0; // k=0 for negative tau (disable update)
    } else if (tau == 0.0) {
      return 1; // k=1 for tau=0
    } else {
      return 1 - Math.exp(-Ts / tau); // fc = 1/(2*pi*tau)
    }
  }

  public static filterCoeffs(
    tau: number,
    Ts: number,
  ): [Tuple<number, 3>, Tuple<number, 2>] {
    const fc = Math.SQRT2 / (2.0 * Math.PI) / tau; // time constant of dampened, non-oscillating part of step response
    const C = Math.tan(Math.PI * fc * Ts);
    const D = C * C + Math.sqrt(2) * C + 1;
    const b0 = (C * C) / D;
    const outB: Tuple<number, 3> = [b0, 2 * b0, b0];
    // a0 = 1.0
    const outA: Tuple<number, 2> = [
      (2 * (C * C - 1)) / D, // a1
      (1 - Math.sqrt(2) * C + C * C) / D, // a2
    ];

    return [outB, outA];
  }

  public static filterInitialState(
    x0: number,
    b: number[],
    a: number[],
  ): Tuple<number, 2> {
    // initial state for steady state (equivalent to scipy.signal.lfilter_zi, obtained by setting y=x=x0 in the filter
    // update equation)
    return [x0 * (1 - b[0]), x0 * (b[2] - a[1])];
  }

  public static filterAdaptStateForCoeffChange<N extends number>(
    last_y: Tuple<number, N>,
    N: N,
    b_old: Tuple<number, 3>,
    a_old: Tuple<number, 2>,
    b_new: Tuple<number, 3>,
    a_new: Tuple<number, 2>,
    state: Tuple<Tuple<number, 2>, N>,
  ): Tuple<Tuple<number, 2>, N> {
    if (Number.isNaN(state[0])) {
      return state;
    }
    for (let i = 0; i < N; i++) {
      state[i][0] = state[i][0] + (b_old[0] - b_new[0]) * last_y[i];
      state[i][1] =
        state[i][1] + (b_old[1] - b_new[1] - a_old[0] + a_new[0]) * last_y[i];
    }

    return state;
  }

  public static filterStep(
    x: number,
    b: Tuple<number, 3>,
    a: Tuple<number, 2>,
    state: Tuple<number, 2>,
  ): [number, Tuple<number, 2>] {
    // difference equations based on scipy.signal.lfilter documentation
    // assumes that a0 == 1.0
    const y = b[0] * x + state[0];
    state[0] = b[1] * x - a[0] * y + state[1];
    state[1] = b[2] * x - a[1] * y;
    return [y, state];
  }

  public static filterVec<N extends number>(
    x: Tuple<number, N>,
    N: N,
    tau: number,
    Ts: number,
    b: Tuple<number, 3>,
    a: Tuple<number, 2>,
    state: Tuple<Tuple<number, 2>, N>,
  ): { state: Tuple<Tuple<number, 2>, N>; out: Tuple<number, N> } {
    // to avoid depending on a single sample, average the first samples (for duration tau)
    // and then use this average to calculate the filter initial state
    const out = tupleOfSize(N, () => 0);
    if (Number.isNaN(state[0])) {
      // initialization phase
      if (Number.isNaN(state[1])) {
        // first sample
        state[0][1] = 0; // state[1] is used to store the sample count
        for (let i = 0; i < N; i++) {
          state[i][0] = 0; // state[2+i] is used to store the sum
        }
      }
      state[0][1]++;
      for (let i = 0; i < N; i++) {
        state[i][0] += x[i];
        out[i] = state[i][0] / state[0][1];
      }
      if (state[0][1] * Ts >= tau) {
        for (let i = 0; i < N; i++) {
          state[i] = VQF.filterInitialState(out[i], b, a);
        }
      }
      return { state, out };
    }

    for (let i = 0; i < N; i++) {
      const [outRes, stateRes] = VQF.filterStep(x[i], b, a, state[i]);
      out[i] = outRes;
      state[i] = stateRes;
    }

    return { state, out };
  }

  public static matrix3SetToScaledIdentity(scale: number): Tuple<number, 9> {
    return [scale, 0, 0, 0, scale, 0, 0, 0, scale];
  }

  public static matrix3Multiply(
    in1: Tuple<number, 9>,
    in2: Tuple<number, 9>,
  ): Tuple<number, 9> {
    return [
      in1[0] * in2[0] + in1[1] * in2[3] + in1[2] * in2[6],
      in1[0] * in2[1] + in1[1] * in2[4] + in1[2] * in2[7],
      in1[0] * in2[2] + in1[1] * in2[5] + in1[2] * in2[8],
      in1[3] * in2[0] + in1[4] * in2[3] + in1[5] * in2[6],
      in1[3] * in2[1] + in1[4] * in2[4] + in1[5] * in2[7],
      in1[3] * in2[2] + in1[4] * in2[5] + in1[5] * in2[8],
      in1[6] * in2[0] + in1[7] * in2[3] + in1[8] * in2[6],
      in1[6] * in2[1] + in1[7] * in2[4] + in1[8] * in2[7],
      in1[6] * in2[2] + in1[7] * in2[5] + in1[8] * in2[8],
    ];
  }

  public static matrix3MultiplyTpsFirst(
    in1: Tuple<number, 9>,
    in2: Tuple<number, 9>,
  ): Tuple<number, 9> {
    return [
      in1[0] * in2[0] + in1[3] * in2[3] + in1[6] * in2[6],
      in1[0] * in2[1] + in1[3] * in2[4] + in1[6] * in2[7],
      in1[0] * in2[2] + in1[3] * in2[5] + in1[6] * in2[8],
      in1[1] * in2[0] + in1[4] * in2[3] + in1[7] * in2[6],
      in1[1] * in2[1] + in1[4] * in2[4] + in1[7] * in2[7],
      in1[1] * in2[2] + in1[4] * in2[5] + in1[7] * in2[8],
      in1[2] * in2[0] + in1[5] * in2[3] + in1[8] * in2[6],
      in1[2] * in2[1] + in1[5] * in2[4] + in1[8] * in2[7],
      in1[2] * in2[2] + in1[5] * in2[5] + in1[8] * in2[8],
    ];
  }

  public static matrix3MultiplyTpsSecond(
    in1: Tuple<number, 9>,
    in2: Tuple<number, 9>,
  ): Tuple<number, 9> {
    return [
      in1[0] * in2[0] + in1[1] * in2[1] + in1[2] * in2[2],
      in1[0] * in2[3] + in1[1] * in2[4] + in1[2] * in2[5],
      in1[0] * in2[6] + in1[1] * in2[7] + in1[2] * in2[8],
      in1[3] * in2[0] + in1[4] * in2[1] + in1[5] * in2[2],
      in1[3] * in2[3] + in1[4] * in2[4] + in1[5] * in2[5],
      in1[3] * in2[6] + in1[4] * in2[7] + in1[5] * in2[8],
      in1[6] * in2[0] + in1[7] * in2[1] + in1[8] * in2[2],
      in1[6] * in2[3] + in1[7] * in2[4] + in1[8] * in2[5],
      in1[6] * in2[6] + in1[7] * in2[7] + in1[8] * in2[8],
    ];
  }

  public static matrix3Inv(input: Tuple<number, 9>): Tuple<number, 9> {
    const A = input[4] * input[8] - input[5] * input[7]; // (e*i - f*h)
    const D = input[2] * input[7] - input[1] * input[8]; // -(b*i - c*h)
    const G = input[1] * input[5] - input[2] * input[4]; // (b*f - c*e)
    const B = input[5] * input[6] - input[3] * input[8]; // -(d*i - f*g)
    const E = input[0] * input[8] - input[2] * input[6]; // (a*i - c*g)
    const H = input[2] * input[3] - input[0] * input[5]; // -(a*f - c*d)
    const C = input[3] * input[7] - input[4] * input[6]; // (d*h - e*g)
    const F = input[1] * input[6] - input[0] * input[7]; // -(a*h - b*g)
    const I = input[0] * input[4] - input[1] * input[3]; // (a*e - b*d)

    const det = input[0] * A + input[1] * B + input[2] * C; // a*A + b*B + c*C;

    if (det >= -EPS && det <= EPS) {
      return tupleOfSize(9, () => 0);
    }

    // out = [A D G; B E H; C F I]/det
    return [
      A / det,
      D / det,
      G / det,
      B / det,
      E / det,
      H / det,
      C / det,
      F / det,
      I / det,
    ];
  }

  protected setup(): void {
    const [outB, outA] = VQF.filterCoeffs(
      this.params.tauAcc,
      this.coeffs.accTs,
    );
    this.coeffs.accLpB = outB;
    this.coeffs.accLpA = outA;

    this.coeffs.kMag = VQF.gainFromTau(this.params.tauMag, this.coeffs.magTs);

    this.coeffs.biasP0 = square(this.params.biasSigmaInit * 100.0);
    // the system noise increases the variance from 0 to (0.1 °/s)^2 in biasForgettingTime seconds
    this.coeffs.biasV =
      (square(0.1 * 100.0) * this.coeffs.accTs) /
      this.params.biasForgettingTime;

    const pMotion = square(this.params.biasSigmaMotion * 100.0);
    this.coeffs.biasMotionW = square(pMotion) / this.coeffs.biasV + pMotion;
    this.coeffs.biasVerticalW =
      this.coeffs.biasMotionW /
      Math.max(this.params.biasVerticalForgettingFactor, 1e-10);

    const pRest = square(this.params.biasSigmaRest * 100.0);
    this.coeffs.biasRestW = square(pRest) / this.coeffs.biasV + pRest;

    const [restGyrB, restGyrA] = VQF.filterCoeffs(
      this.params.restFilterTau,
      this.coeffs.gyrTs,
    );
    this.coeffs.restGyrLpB = restGyrB;
    this.coeffs.restGyrLpA = restGyrA;
    const [restAccB, restAccA] = VQF.filterCoeffs(
      this.params.restFilterTau,
      this.coeffs.accTs,
    );
    this.coeffs.restAccLpB = restAccB;
    this.coeffs.restAccLpA = restAccA;

    this.coeffs.kMagRef = VQF.gainFromTau(
      this.params.magRefTau,
      this.coeffs.magTs,
    );
    if (this.params.magCurrentTau > 0) {
      const [magB, magA] = VQF.filterCoeffs(
        this.params.magCurrentTau,
        this.coeffs.magTs,
      );
      this.coeffs.magNormDipLpB = magB;
      this.coeffs.magNormDipLpA = magA;
    } else {
      this.coeffs.magNormDipLpB.fill(NaN);
      this.coeffs.magNormDipLpA.fill(NaN);
    }

    this.resetState();
  }

  protected params: VQFParams = defaultParams();

  protected state: VQFState = defaultState();

  protected coeffs: VQFCoefficients = defaultCoeffs();
}
