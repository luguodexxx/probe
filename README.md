## ProbeDesign: A RNA probe design tools [![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

<pre>
__________              ___.          ________                .__               
\______   \_______  ____\_ |__   ____ \______ \   ____   _____|__| ____   ____  
 |     ___/\_  __ \/  _ \| __ \_/ __ \ |    |  \_/ __ \ /  ___/  |/ ___\ /    \ 
 |    |     |  | \(  <_> ) \_\ \  ___/ |    `   \  ___/ \___ \|  / /_/  >   |  \
 |____|     |__|   \____/|___  /\___  >_______  /\___  >____  >__\___  /|___|  /
                             \/     \/        \/     \/     \/  /_____/      \/ 
</pre>


ProbeDesign, a simple command line to generate probe. It's  ~~still not~~  completed.


### INSTALL
python >= 3.6
```shell
python setup.py install
```
### UASGE
There are two function.

#### transcript

Generate the probe based on the transcript ID.

```shell
probedesign transcripts -index TranscriptIndex -ta example/test.id
```

#### junction

Generate the probe based on the junction information.

```shell
probedesign junction -faG GenomeFastaFile -ta example/test.junction -index GenomeIndex.bowtie2
```

#### circ

Generate the probe based on the circ-junction information.

```shell
probedesign circ -faG GenomeFastaFile -ta example/test.circRNA -index GenomeIndex.bowtie2
```


### Reference
We generate the probe followed the same filter conditions with blockParse.py in [OligoMiner](https://github.com/brianbeliveau/OligoMiner)
