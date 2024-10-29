function result = jelena_ecg_helper(x, qrs_i_raw, ECG)
    temp = qrs_i_raw + x;
    temp = temp(temp > 0 & temp < size(ECG.data,2));
    result = mean(ECG.data(1,temp));
end
