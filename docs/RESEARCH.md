# Research on Solid-State Transformers (SST) and Silicon Carbide (SiC)

This document summarizes updated research for the Compact In-Line Solid-State Substation (CSSS) project, specifically focusing on the implementation of high-frequency solid-state transformers using Silicon Carbide (SiC) technology.

## High-Voltage SiC MOSFETs
Recent advancements by companies like Wolfspeed have made 10kV SiC MOSFETs commercially viable. These high-voltage devices are crucial for direct AC-AC conversion from transmission levels (e.g., 138kV) to distribution levels. The 10kV SiC MOSFETs enable simplified topologies with fewer series-connected cells compared to conventional silicon IGBT-based designs. Furthermore, they exhibit excellent bipolar stability, addressing a major failure mode in high-voltage SiC devices over long periods of operation.

## High-Frequency Operation
A key benefit of SiC technology is its ability to operate efficiently at high switching frequencies. For compact magnetics, SST power stages typically target switching frequencies higher than 10 kHz. Operating above 10 kHz directly reduces the size and weight of the magnetic components (like the autotransformer core in the CSSS design), because the required magnetic coils and core size decrease as frequency increases. Wolfspeed's 10kV SiC MOSFETs can achieve this with substantially lower switching loss and faster transitions compared to conventional silicon devices, making compact, lightweight SSTs feasible.

## Core Characteristics and Implications
Operating an amorphous magnetic core at these elevated frequencies (>10 kHz) requires careful consideration of non-linear B-H curve behavior to manage saturation and core loss calculations. The combination of high-frequency operation and a direct AC-AC topology (such as a Matrix Converter or AC-AC Dual Active Bridge) minimizes components and conversion losses, fulfilling the "extreme decisions" for maximum compactness in the CSSS architecture.
