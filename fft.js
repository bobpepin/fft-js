export function dft(inputRe, inputIm, outputRe, outputIm) {
    const N = inputRe.length;
    for(let k=0; k < N; k++) {
        outputRe[k] = outputIm[k] = 0;
        for(let n=0; n < N; n++) {
            const twiddleRe = Math.cos(-2*Math.PI*k*n/N);
            const twiddleIm = Math.sin(-2*Math.PI*k*n/N);
            outputRe[k] += twiddleRe*inputRe[n] - twiddleIm*inputIm[n];
            outputIm[k] += twiddleIm*inputRe[n] + twiddleRe*inputIm[n];
        }
    }
}

function dftStride(inputRe, inputIm, outputRe, outputIm, N,
             inputOffset, inputStride, outputOffset, outputStride) {
    for(let k=0; k < N; k++) {
        const ki = outputOffset + outputStride*k;        
        outputRe[ki] = outputIm[ki] = 0;
        for(let n=0; n < N; n++) {
            const ni = inputOffset + inputStride*n;
            const twiddleRe = Math.cos(-2*Math.PI*k*n/N);
            const twiddleIm = Math.sin(-2*Math.PI*k*n/N);
            outputRe[ki] += twiddleRe*inputRe[ni] - twiddleIm*inputIm[ni];
            outputIm[ki] += twiddleIm*inputRe[ni] + twiddleRe*inputIm[ni];
        }
    }
}

export function dft2d(dimensions, inputRe, inputIm, buffersRe, buffersIm) {
    const [dim1, dim2] = dimensions;
    for(let i=0; i < dim1; i++) {
        let offset = i*dim2;
        let stride = 1;
        dftStride(inputRe, inputIm, buffersRe[0], buffersIm[0], dim2,
                  offset, stride, offset, stride);
    }
    for(let i=0; i < dim2; i++) {
        const offset = i;
        const stride = dim1;
        dftStride(
            buffersRe[0], buffersIm[0], buffersRe[1], buffersIm[1],
            dim1, offset, stride, offset, stride
        );
    }
    return 1;
}

export function fft1(inputRe, inputIm,  
                     N, stride, ofs, 
                     twiddleRe, twiddleIm, twiddleStride,
                     fftOffset,
                     outputStride, outputOffset, outputRe, outputIm) {
    if(N == 1) {
        outputRe[fftOffset*outputStride + outputOffset] = inputRe[ofs];
        outputIm[fftOffset*outputStride + outputOffset] = inputIm[ofs];
    } else {
//         dftStride(inputRe, inputIm, outputRe, outputIm, N/2,
//              ofs, 2*stride, outputOffset, 1);
//         dftStride(inputRe, inputIm, outputRe, outputIm, N/2,
//              stride+ofs, 2*stride, outputOffset+N/2, 1);
//         console.log(`[${outputRe}]`, `[${outputIm}]`);
        fft1(inputRe, inputIm, N/2, 2*stride, ofs, twiddleRe, twiddleIm, twiddleStride, fftOffset, outputStride, outputOffset, outputRe, outputIm);
        fft1(inputRe, inputIm, N/2, 2*stride, stride+ofs, twiddleRe, twiddleIm, twiddleStride, fftOffset + N/2, outputStride, outputOffset, outputRe, outputIm);
        for(let k=0; k < N/2; k++) {
            const evenIndex = (fftOffset+k)*outputStride + outputOffset;
            const oddIndex = (fftOffset+N/2+k)*outputStride + outputOffset;
            const outputIndex1 = evenIndex;
            const outputIndex2 = oddIndex;
            const twiddleIndex = k*twiddleStride + (N-1);
            const twiddledOddRe = twiddleRe[twiddleIndex]*outputRe[oddIndex] - twiddleIm[twiddleIndex]*outputIm[oddIndex];
            const twiddledOddIm = twiddleIm[twiddleIndex]*outputRe[oddIndex] + twiddleRe[twiddleIndex]*outputIm[oddIndex];
//             console.log(`Before: [${outputRe}]`);
            const evenRe = outputRe[evenIndex];
            const evenIm = outputIm[evenIndex];
//             console.log("k", k, "N", N, 
//                         "twiddle val", twiddleRe[twiddleIndex], twiddleIm[twiddleIndex], 
//                         "even val", evenRe, evenIm,
//                         "odd val", outputRe[oddIndex], outputIm[oddIndex],
//                         "twiddled", twiddledOddRe, twiddledOddIm);            
            outputRe[outputIndex1] = evenRe + twiddledOddRe;
            outputIm[outputIndex1] = evenIm + twiddledOddIm;
            outputRe[outputIndex2] = evenRe - twiddledOddRe;
            outputIm[outputIndex2] = evenIm - twiddledOddIm;
//             console.log(`After: [${outputRe}]`);            
        }
    }
}

export function computeTwiddleFactors(Nmax, invert, outputRe, outputIm) {
    for(let k=0; k < Nmax; k++) {
        for(let N=1; N <= Nmax; N++) {
            outputRe[k*Nmax+(N-1)] = Math.cos(-2*invert*Math.PI*k/N);
            outputIm[k*Nmax+(N-1)] = Math.sin(-2*invert*Math.PI*k/N);
        }
    }
}

export function fft1d(inputRe, inputIm, outputRe, outputIm, twiddleRe, twiddleIm, twiddleStride) {
    const N = inputRe.length;
    if(twiddleRe === undefined) {
        twiddleRe = new Float32Array(N*N);
        twiddleIm = new Float32Array(N*N);
        computeTwiddleFactors(N, 1, twiddleRe, twiddleIm);
        twiddleStride = N;
    }
    console.log("twiddle", twiddleRe, twiddleIm);
    fft1(inputRe, inputIm, N, 1, 0, twiddleRe, twiddleIm, twiddleStride, 0, 1, 0, outputRe, outputIm);
}

// export function fft1(inputRe, inputIm,  
//                      N, stride, ofs, 
//                      twiddleRe, twiddleIm, twiddleStride,
//                      fftOffset
//                      outputStride, outputOffset, outputRe, outputIm) {



export function fft2d(dimensions, inputRe, inputIm, buffersRe, buffersIm, twiddleRe, twiddleIm, twiddleStride) {
    const [dim1, dim2] = dimensions;
    if(twiddleRe === undefined) {
        twiddleStride = Math.max(...dimensions);
        twiddleRe = new Float32Array(twiddleStride*twiddleStride);
        twiddleIm = new Float32Array(twiddleStride*twiddleStride);
        computeTwiddleFactors(twiddleStride, 1, twiddleRe, twiddleIm);
    }
    for(let i=0; i < dim1; i++) {
        let offset = i*dim2;
        let stride = 1;
        fft1(inputRe, inputIm, dim2, stride, offset, twiddleRe, twiddleIm, twiddleStride, 0, stride, offset, buffersRe[0], buffersIm[0]);
    }
    for(let i=0; i < dim2; i++) {
        const offset = i;
        const stride = dim1;
        fft1(
            buffersRe[0], buffersIm[0], dim1, stride, offset, twiddleRe, twiddleIm, twiddleStride, 0, stride, offset, buffersRe[1], buffersIm[1]
        );
    }
    return 1;
}

