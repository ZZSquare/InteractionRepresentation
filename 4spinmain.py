import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from qiskit.providers.aer import noise
from qiskit.ignis.mitigation.measurement import complete_meas_cal, CompleteMeasFitter
from qiskit.transpiler import Layout
from qiskit import Aer, QuantumRegister, QuantumCircuit, execute, IBMQ, ClassicalRegister
from qfunctions import adiabaticramp, theta, fourspin_instruction, exactm, interactionm


def twospin_df(J, Bx, dB, Bz_max, dt_steps, dt_steps_bool, gamma, skip, provider, backend, pq1, pq2):
    count_list = []
    calibrated_count_list = []
    tarray, Bzarray, dt_or_steps, _ = adiabaticramp(J, Bx, dB, Bz_max, dt_steps, dt_steps_bool, gamma)

    if dt_steps_bool == 'dt':
        dt = dt_steps
    else:
        dt = dt_or_steps

    thetaarray, _ = theta(Bzarray, dt)

    # Calibration Matrix
    meas_calibs, state_labels = complete_meas_cal(qubit_list=[0, 1, 2, 3], qr=QuantumRegister(4), circlabel='4spin')
    cal_results = execute(meas_calibs, backend=Aer.get_backend('qasm_simulator'), shots=1000).result()
    meas_calibs, state_labels = complete_meas_cal(qubit_list=[0, 1, 2, 3], qr=QuantumRegister(4), circlabel='4spin')
    meas_fitter = CompleteMeasFitter(cal_results, state_labels, circlabel='4spin')
    #meas_fitter.plot_calibration()

    qmap = Layout()

    i = 0
    while (i < len(thetaarray)):
        print('Bz = %f' % (Bzarray[i]))
        qr = QuantumRegister(4, 'qr')
        cr = ClassicalRegister(4, 'cr')
        circ = QuantumCircuit(qr, cr)
        #qmap.from_dict(input_dict={qr[0]: pq1, qr[1]: pq2})
        circ.initialize([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 1, 2, 3])

        fourspin_instruction(circ, J, Bx, thetaarray[:i + 1], dt)
        circ.measure([0, 1, 2, 3], [0, 1, 2, 3])

        result = execute(circ, Aer.get_backend('qasm_simulator'), shots=5000).result()
        counts = result.get_counts()
        counts['Time'] = tarray[i]
        counts['Bz'] = Bzarray[i]
        count_list.append(counts)
        with open("count_list.txt", "w") as output:
            output.write(str(count_list))

        mitigated_counts = meas_fitter.filter.apply(result).get_counts()
        mitigated_counts['Time'] = tarray[i]
        mitigated_counts['Bz'] = Bzarray[i]
        calibrated_count_list.append(mitigated_counts)
        with open("calibrated_count_list.txt", "w") as output:
            output.write(str(calibrated_count_list))

        i = i + 1 + skip

    # Creating dataframe
    df = pd.DataFrame(count_list)
    time_col = df.pop('Time')
    df.insert(0, 'Time', time_col)
    df['Exact'] = exactm(J, Bx, Bzarray, dt)
    df['Interaction'] = interactionm(J, Bx, thetaarray, dt)

    calibrated_df = pd.DataFrame(calibrated_count_list)
    time_col = df.pop('Time')
    df.insert(0, 'Time', time_col)

    if '0000' not in df:
        df['0000'] = 0
    if '0001' not in df:
        df['0001'] = 0
    if '0010' not in df:
        df['0010'] = 0
    if '0011' not in df:
        df['0011'] = 0
    if '0100' not in df:
        df['0000'] = 0
    if '0101' not in df:
        df['0001'] = 0
    if '0110' not in df:
        df['0010'] = 0
    if '0111' not in df:
        df['0011'] = 0
    if '1000' not in df:
        df['0000'] = 0
    if '1001' not in df:
        df['0001'] = 0
    if '1010' not in df:
        df['0010'] = 0
    if '1011' not in df:
        df['0011'] = 0
    if '1100' not in df:
        df['0000'] = 0
    if '1101' not in df:
        df['0001'] = 0
    if '1110' not in df:
        df['0010'] = 0
    if '1111' not in df:
        df['0011'] = 0
    df = df.fillna(0)

    if '0000' not in calibrated_df:
        calibrated_df['00'] = 0
    if '0001' not in calibrated_df:
        calibrated_df['01'] = 0
    if '0010' not in calibrated_df:
        calibrated_df['10'] = 0
    if '0011' not in calibrated_df:
        calibrated_df['11'] = 0
    if '0100' not in calibrated_df:
        calibrated_df['00'] = 0
    if '0101' not in calibrated_df:
        calibrated_df['01'] = 0
    if '0110' not in calibrated_df:
        calibrated_df['10'] = 0
    if '0111' not in calibrated_df:
        calibrated_df['11'] = 0
    if '1000' not in calibrated_df:
        calibrated_df['00'] = 0
    if '1001' not in calibrated_df:
        calibrated_df['01'] = 0
    if '1010' not in calibrated_df:
        calibrated_df['10'] = 0
    if '1011' not in calibrated_df:
        calibrated_df['11'] = 0
    if '1100' not in calibrated_df:
        calibrated_df['00'] = 0
    if '1101' not in calibrated_df:
        calibrated_df['01'] = 0
    if '1110' not in calibrated_df:
        calibrated_df['10'] = 0
    if '1111' not in calibrated_df:
        calibrated_df['11'] = 0
    calibrated_df = calibrated_df.fillna(0)

    # Calculating mz
    total = df['0000'] + df['0001'] + df['0010'] + df['0011'] + df['0100'] + df['0101'] + df['0110'] + df['0111'] + df[
        '1000'] + df['1001'] + df['1010'] + df['1011'] + df['1100'] + df['1101'] + df['1110'] + df['1111']
    df['mz'] = -(2*df['0000'] + df['0001'] + df['0010'] + df['0100'] + df['1000'] - df['0111'] - df['1011'] - df['1101'] - df['1110'] - 2*df['1111']) / total
    calibrated_df['mz'] = -(2*df['0000'] + df['0001'] + df['0010'] + df['0100'] + df['1000'] - df['0111'] - df['1011'] - df['1101'] - df['1110'] - 2*df['1111']) / total

    # Creating Files
    if dt_steps_bool == 'dt':
        df.to_csv('4spin J={:.3f} Bx={:.3f} dB={:.3f} BzMax={:.3f} dt={:.3f} gamma={:.3f}.csv'.format(J, Bx, dB, Bz_max,
                                                                                                      dt_steps, gamma))
        calibrated_df.to_csv(
            'Calibrated 4spin J={:.3f} Bx={:.3f} dB={:.3f} BzMax={:.3f} dt={:.3f} gamma={:.3f}.csv'.format(J, Bx, dB,
                                                                                                           Bz_max,
                                                                                                           dt_steps,
                                                                                                           gamma))
    else:
        df.to_csv(
            '4spin J={:.3f} Bx={:.3f} dB={:.3f} BzMax={:.3f} BzSteps={:.3f} gamma={:.3f}.csv'.format(J, Bx, dB, Bz_max,
                                                                                                     dt_steps, gamma))
        calibrated_df.to_csv(
            'Calibrated 4spin J={:.3f} Bx={:.3f} dB={:.3f} BzMax={:.3f} BzSteps={:.3f} gamma={:.3f}.csv'.format(J, Bx,
                                                                                                                dB,
                                                                                                                Bz_max,
                                                                                                                dt_steps,
                                                                                                                gamma))

    return df, calibrated_df, dt_or_steps, thetaarray

def showfig(J, Bx, dB, Bz_max, dt_steps, dt, gamma, skip, provider, backend, pq1, pq2):
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    df, calibrated_df, _, thetaarray = twospin_df(J, Bx, dB, Bz_max, dt_steps, dt, gamma, skip, provider, backend, pq1, pq2)

    fig.add_trace(go.Scatter(x=calibrated_df['Bz'], y=df['mz'], name='Calibrated', mode='lines'), secondary_y=False)
    fig.add_trace(go.Scatter(x=df['Bz'], y=df['mz'], name='Circuit', mode='lines'), secondary_y=False)
    fig.add_trace(go.Scatter(x=df['Bz'], y=df['Exact'], name='Exact', mode='lines'), secondary_y=False)
    fig.add_trace(go.Scatter(x=df['Bz'], y=df['Interaction'], name='Interaction', mode='lines'), secondary_y=False)

    # fig.add_trace(go.Scatter(x=df['Bz'], y=thetaarray, name='Theta', mode='markers'), secondary_y=True)
    title = "J=" + str(J) + " Bx=" + str(Bx) + " gamma=" + str(gamma) + " dt_steps=" + str(dt_steps)
    fig.update_layout(title_text=title)
    fig.update_xaxes(title_text="Bz")
    fig.update_yaxes(title_text="mz", secondary_y=False)
    # fig.update_yaxes(title_text="theta", secondary_y=True)
    fig.show()



def main():

    # Set account credentials here
    #IBMQ_KEY =
    #IBMQ.save_account(IBMQ_KEY)
    provider = IBMQ.load_account()

    # Set name of backend
    # eg. 'ibmq_santiago'
    backend ='ibmq_santiago'

    # Map virtual qubits to physical qubits
    physicalq1 =0
    physicalq2 =7

    J = 1
    Bx = 0.05
    dB = 0.01
    Bz_max = 2
    dt_steps = 0.5
    dt = 'dt'
    gamma = 10

    showfig(J, Bx, dB, Bz_max, dt_steps, dt, gamma, 0, provider, backend, physicalq1, physicalq2)

if __name__ == "__main__":
    main()