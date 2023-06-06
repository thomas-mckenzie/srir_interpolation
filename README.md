# Interpolation of Spatial Room Impulse Responses for Room Transitions

This repository contains MATLAB functions to interpolate between spatial room impulse responses at different receiver positions, suitable for room transitions and other non-simple recording setups. 
Please [refer to the publication](https://www.aes.org/e-lib/browse.cfm?elib=22140) for more information:
   ```
T. McKenzie, N. Meyer-Kahlen, C. Hold, S. J. Schlecht, and V. Pulkki: "Auralisation of measured room transitions in virtual
reality." Journal of the Audio Engineering Society, doi: 10.17743/jaes.2022.0084, 2023.
   ```  
   ```
T. McKenzie, N. Meyer-Kahlen, R. Daugintis, L. McCormack, S. J. Schlecht, and V. Pulkki: "Perceptually informed interpolation and  
rendering of spatial room impulse responses for room transitions." International Congress on Acoustics, Gyeongju, pp. 1-11. 2022.
   ```

   
Two demo scripts 'demo_interpolate_SRIRs_function' are included: one for 1D (ie along a line) measurements for the transition between two coupled rooms, and one for 2D measurements inside a single room. The script also includes an example of formatting for saving a sofa file of the resulting SRIRs. 

For listening and auralising sofa files in this format, please see the [Sparta 6DOFConv VST plugin](https://leomccormack.github.io/sparta-site/docs/plugins/sparta-suite/#6dofconv). 

# Source Interpolation of Spatial Room Impulse Responses

This repository also contains MATLAB functions to interpolate between spatial room impulse responses at different source positions. 
Please [refer to the publication](https://www.aes.org/e-lib/browse.cfm?elib=22044) for more information:
   
   ```
T. McKenzie, S. J. Schlecht: "Source position interpolation of spatial room impulse responses" AES 154th Convention, Espoo, 
Helsinki, Finland, 2023. 
   ```
   
A demonstration script 'demo_interpolate_SRIRs_source.m' is included, which along with 'analyse_interpolate_SRIRs_source_AES.m' plots the figures from the AES convention paper. An example of formatting for saving a sofa file of the resulting SRIRs is included, which is compatible with the [Sparta 6DoFConv VST plugin](https://leomccormack.github.io/sparta-site/docs/plugins/sparta-suite/#6dofconv), for listening and auralising up to 6DoF sofa files.
