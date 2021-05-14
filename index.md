---
layout: workshop
venue: "University at Buffalo, SUNY"   # brief name of host site without address (e.g., "Euphoric State University")
address: "University at Buffalo, SUNY, North Campus"     # full street address of workshop (e.g., "Room A, 123 Forth Street, Blimingen, Euphoria")
country: "United States"      # lowercase two-letter ISO country code such as "fr" (see https://en.wikipedia.org/wiki/ISO_3166-1#Current_codes)
language: "en"     # lowercase two-letter ISO language code such as "fr" (see https://en.wikipedia.org/wiki/List_of_ISO_639-1_codes)
latitude: "43.002890"     # decimal latitude of workshop venue (use https://www.latlong.net/)
longitude: "-78.788780"    # decimal longitude of the workshop venue (use https://www.latlong.net)
humandate: "Jun 14-26, 2021"    # human-readable dates for the workshop (e.g., "Feb 17-18, 2020")
humantime: "11:00 am - 5:00 pm EDT"    # human-readable times for the workshop (e.g., "9:00 am - 4:30 pm")
startdate: 2021-06-14      # machine-readable start date for the workshop in YYYY-MM-DD format like 2015-01-01
enddate: 2021-06-26        # machine-readable end date for the workshop in YYYY-MM-DD format like 2015-01-02
instructor: ["Alexey Akimov", "Jeanette Sperhac", "Ivan Infante",
             "Felipe Zapata", "Sergei Tretiak", "Walter Malone", 
             "Mario Barbatti", "Hans Lischka", "Amber Jain"]  # boxed, comma-separated list of instructors' names as strings, like ["Kay McNulty", "Betty Jennings", "Betty Snyder"]
helper: ["Mohammad Shakiba"]     # boxed, comma-separated list of helpers' names, like ["Marlyn Wescoff", "Fran Bilas", "Ruth Lichterman"]
email: ["alexeyak@buffalo.edu"]    # boxed, comma-separated list of contact email addresses for the host, lead instructor, or whoever else is handling questions, like ["marlyn.wescoff@example.org", "fran.bilas@example.org", "ruth.lichterman@example.org"]
collaborative_notes:         # optional: URL for the workshop collaborative notes, e.g. an Etherpad or Google Docs document (e.g., https://pad.carpentries.org/2015-01-01-euphoria)
googleform: https://forms.gle/kup1bkWibLsHH8Bn7
carpentry: "sc"
---


{% comment %}
root: .  # Is the only page that doesn't follow the pattern /:path/index.html
permalink: index.html  # Is the only page that doesn't follow the pattern /:path/index.html
{% endcomment %}


> # UPDATE: The CyberTraining 2021 summer school and workshop registration is now open.
> Note, the workshop is offered in an online format.
> {: .testimonial}

# CyberTraining 2021: Excited States and Nonadiabatic Dynamics Summer School 

## About the Summer School and Workshop

The CyberTraining workshop aims to educate graduate students, postdocs, researchers, and educators working in a 
broader field of nonadiabatic and excited-state dynamics as well as in computational material sciences in a variety
of tools and methods for such types of calculations. The workshop will provide conceptual and practical hands-on 
training in a range of methods and cyberinfrastructure (software and platforms) for modeling excited state and 
nonadiabatic dynamics in abstract models and atomistic materials. We will also cover tools and workflows for 
building atomistic models, computing excited states of molecular and periodic systems, as well as pre- and post-processing 
operations, and data analysis. 

Participants will not only learn about using the tools but will be exposed to the underlying machinery 
of such methods and will be familiarized with their development. The programming-driven nature of the school will
help the participants to go beyond the standard computational chemistry curriculum. The workshop will culminate with 
a capstone project presentation, through which the participants will demonstrate their ability to leverage the new tools 
in their active research. 

Keywords and topics:

- nonadiabatic dynamics
- excited states
- quantum dynamics
- quantum-classical methods
- charge transfer
- excitation energy transfer
- trajectory surface hopping
- TD-DFT
- algorithms and methods
- software, programming, Python
- best practices, git, github

 
The school aims to provide training in a range of advanced tools. 
This year, the focus will be on the following packages:

- Libra (Akimov)
- NEXMD (Tretiak)
- Newton-X (Barbatti)
- nano-qmflows (Infante, Zapata)
- CAT, auto-FOX (Infante, Zapata)
- COLUMBUS (Lischka)
- DFTB+
- CP2K
- Quantum Espresso
- ErgoSCF

The school will leverage the OnDemand gateway at the University at Buffalo


## Logistics

{% if page.humandate %}
<p id="when">
  <strong>When:</strong>
  {{page.humandate}}.
  {% include workshop_calendar.html %}
</p>
{% endif %}

{% if page.latitude and page.longitude %}
<p id="where">
  <strong>Where:</strong>
  {{page.address}}.
  Get directions with
  <a href="//www.openstreetmap.org/?mlat={{page.latitude}}&mlon={{page.longitude}}&zoom=16">OpenStreetMap</a>
  or
  <a href="//maps.google.com/maps?q={{page.latitude}},{{page.longitude}}">Google Maps</a>.
</p>
{% endif %}

{% comment %}
CONTACT EMAIL ADDRESS
Display the contact email address set in the configuration file.
{% endcomment %}
<p id="contact">
  <strong>Contact</strong>:
  Please email
  {% if page.email %}
  {% for email in page.email %}
  {% if forloop.last and page.email.size > 1 %}
  or
  {% else %}
  {% unless forloop.first %}
  ,
  {% endunless %}
  {% endif %}
  <a href='mailto:{{email}}'>{{email}}</a>
  {% endfor %}
  {% else %}
  to-be-announced
  {% endif %}
  for more information.
</p>

### Schedule

{% include base_path.html %}

The details may vary and the order of topics may be changed, the topics may be omitted or added. Please check for the updates. 

  <table class="table table-striped">
  <tr>
    <td class="col-md-3">June 14, 2021 (Day 1), <strong>Monday</strong></td>
    <td class="col-md-7">
    * Introduction. Overview of CyberInfrastructure. 
    * Revision of Python and C++ programming. Coding Molecular Dynamics and Monte Carlo.
    </td>
    <td class="col-md-2">Schedule of Day 1</td>
  </tr>

  <tr>
    <td class="col-md-3">August 3, 2020 (Day 2), Monday</td>
    <td class="col-md-7">Opening. <a href="01-introduction">Introduction. VIDIA gateway. Setting up.  
    Revision of Python programming. Revision of the best practices</a></td>
    <td class="col-md-2">Schedule of Day 2</td>      
  </tr>

  <tr>
    <td class="col-md-3">August 4, 2020 (Day 3), Tuesday</td>
    <td class="col-md-7"><a href="02-methodology-libra">Methodology prototyping in Python. 
    Introduction of the Libra package. Classical molecular dynamics and convenience functions.</a></td>
    <td class="col-md-2">Schedule of Day 3</td>
  </tr>


  <tr>
    <td class="col-md-3">August 5, 2020 (Day 4), Wednesday</td>
    <td class="col-md-7">Wavepacket dynamics on a grid. Hierarchy of equations of motion (HEOM).</td>
    <td class="col-md-2">Schedule of Day 4</td>
  </tr>

  <tr>
    <td class="col-md-3">August 6, 2020 (Day 5), Thursday</td>
    <td class="col-md-7">Trajectory surface hopping and Ehrenfest methods for model Hamiltonians.</td>
    <td class="col-md-2">Schedule of Day 5</td>
  </tr>

  <tr>
    <td class="col-md-3">August 7, 2020 (Day 6), Friday</td>
    <td class="col-md-7">Nonadiabatic dynamics of atomistics systems: Libra/ErgoSCF, Libra/DFTB+, and Libra/QuantumEspresso.</td>
    <td class="col-md-2">Schedule of Day 6</td>
  </tr>

  <tr>
    <td class="col-md-3">August 8, 2020 (Day 7), Saturday</td>
    <td class="col-md-7">On your own. Projects time</td>
    <td class="col-md-2">Schedule of Day 7</td>
  </tr>

  <tr>
    <td class="col-md-3">August 9, 2020 (Day 8), <strong>Sunday</strong></td>
    <td class="col-md-7">On your own. Projects time</td>
    <td class="col-md-2">Schedule of Day 8</td>
  </tr>


  <tr>
    <td class="col-md-3">August 10, 2020 (Day 9), Monday</td>
    <td class="col-md-7">NA-MD of nanoclusters: Hands on with QMflows-NAMD and CP2K</td>
    <td class="col-md-2">Schedule of Day 9</td>
  </tr>

  <tr>
    <td class="col-md-3">August 11, 2020 (Day 10), Tuesday</td>
    <td class="col-md-7">NA-MD of nanoclusters and condensed matter: Hands on with Libra workflows, DFTB+ and Quantum Espresso</td>
    <td class="col-md-2">Schedule of Day 10</td>
  </tr>

  <tr>
    <td class="col-md-3">August 12, 2020 (Day 11), Wednesday</td>
    <td class="col-md-7">NA-MD of large organic molecules: Hands on with NEXMD</td>
    <td class="col-md-2">Schedule of Day 11</td>
  </tr>

  <tr>
    <td class="col-md-3">August 13, 2020 (Day 12), Thursday</td>
    <td class="col-md-7">Wavepacket dynamics tutorial. NA-MD of small molecules: Hands on with SHARC and OpenMolcas</td>
    <td class="col-md-2">Schedule of Day 12</td>
  </tr>

  <tr>
    <td class="col-md-3">August 14, 2020 (Day 13), Friday</td>
    <td class="col-md-7"> ... Continued </td>
    <td class="col-md-2">Schedule of Day 13</td>
  </tr>

  <tr>
    <td class="col-md-3">August 15, 2020 (Day 14), Saturday</td>
    <td class="col-md-7">Project presentations</td>
    <td class="col-md-2">Schedule of Day 14</td>
  </tr>

  <tr>
    <td class="col-md-3">August 16, 2020 (Day 15), <strong>Sunday</strong></td>
    <td class="col-md-7">Depart</td>
    <td class="col-md-2"></td>
  </tr>
  </table>



## Participation
### How to apply to the school

1. Read this page carefully
2. Prepare your application package (you will need it in the next steps)

   2.1. your CV (including graduate GPA)

   2.2. a statement of purpose PDF should describe in no more than 2 pages:

   * your current/ongoing research projects and interests; 
   * how you plan to use the CyberTraining skill gained in this workshop in your research, for instance if you expect using any of the
     packages that will be covered at this workshop (see the agenda);
   * propose at least one potential project to be completed during the summer school; the project will be presented at the end of the 
     event and should involving one or more tools/software covered during the workshop (see the agenda). The quality and feasibility 
     of the proposed workshop projects will be considered during the selection of the participants. 
         
   2.3. request your advisor to submit a letter of recommendation for you to the following email: "alexeyak AT buffalo DOT edu", 
   please replace "AT" and "DOT" with the corresponding characters

3. Complete the <a href="https://forms.gle/kup1bkWibLsHH8Bn7" target="_blank" rel="nofollow">**Registration form**</a>


### Important dates
   * School application materials are due 5 pm EDT, June 7, 2020
   * Students and Postdocs will be notified of their admission by June 10, 2020
   * School starts: 11 am EDT, June 14, 2020


### Who can apply

This summer school is primarily for graduate students working in computational 
modeling of excited states and nonadiabatic dynamics, both in abstract and atomistic
applications/problems. The school aims to help researchers/students working either in 
methodology development for nonadiabatic or quantum-classical dynamics and in 
applied studies of various types of solar energy materials (photovoltaics, photocatalytics, etc.). 

Postdocs and researchers wishing to acquire the practical experience with new simulation
tools and expand their knowledge in the areas of excited states and nonadiabatic dynamics
are also welcomed to participate.




### Selection and restrictions

* **Competitive selection** The applicants will be selected based on the strength of their statement of purpose, as well as the adequate 
  support of their supervisors and their level of fundamental preparation. The lack of training in specialized methods and software is not a problem. 
  What is more important is how ready the applicants are to absorb the new knowlege, how efficiently they can operate during the workshop, 
  and how critical the use of the methods/tools covered in the workshop may be for your future research or career (e.g. educating others). 

* **The VPN cap/instructing efficiency limit.** The hands-on session will be facilitated by the CyberInfrastructure built
  at the UB CCR cluster. As such, users have to use VPN to remotely access the cluster. The UBIT department has provided 
  a block of 30 external VPNs for non-UB participants (including about 10 instructors). This number sets the limit of about 
  20 people for non-instructor participants we can accept to the fully-fledged (talks/demos + hands-on) event. However, more 
  participants may be admitted to the theory talks/demo sessions. 

* **Export control limit.** Certain countries (e.g. China, Iran, Russia, etc.) can not be issued the UB VPN, so the participants 
  from these countries can not use the UB CCR cyberinftrastructure during the hands-on activities. Such participants may still be 
  admitted to the theory talks/demo sessions. 

* **Group champions.** We anticipate the workshop may be of interest to more than 1 person from any given research group. To broaden 
  and diversify the participation, we will admit only 1 person from any research group (2, if we have room) to a fully-fledged workshop
  (talks/demos + hands-on). This is the group champion. Although this person will not be allowed to share their login credentials with other
  group members, they are free to communicate with other group members (that may be accepted to the workshop as non-champions) during hands-on
  exercises and share their screens with the group mates.





### Acknowledgement

This workshop is made possible by the NSF-OAC CyberTraining program. Thank you!



{% include links.md %}
